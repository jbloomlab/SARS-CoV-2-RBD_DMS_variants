"""Take Excel file from plate reader and conver to fraction infectivity."""


import argparse
import itertools
import os

import numpy as np

import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Convert plate reader '
                                     'Excel file to fraction infectivity '
                                     'csv.')
    parser.add_argument('infile', type=str, help='Path to excel file '
                        'to convert to fraction infectivity.')
    parser.add_argument('outfile', type=str, help="Path for output "
                        "fraction infectivity csvs.")
    parser.add_argument('sheet_map', type=str, help="File to map "
                        "plate number and samples. Must have columns:"
                        "'Plate', 'Sample', 'Virus', 'SampleNum', "
                        "'PlateLayout', 'StartDil', and 'DilFactor'")
    parser.add_argument('plate_layouts_dir', type=str, help='Directory '
                        'containing csv files specifying plate layouts.')
    parser.add_argument('--allow_cells_only_bg',
                        help="If no 'neg_ctrl', do we use 'cells_only' for "
                             "quantify background?",
                        action='store_true',
                        )
    return parser.parse_args()


def get_locs(layout, value):
    """Get (index, column) tuples for location of value in layout df."""
    locs_list = []
    locs = layout.isin([value])

    series = locs.any()
    columns = list(series[series].index)

    for col in columns:
        rows = list(locs[col][locs[col]].index)

        for row in rows:
            locs_list.append((row, col))

    return locs_list


def main():
    """Write fraction infectivity csv from plate reader data.

    Use user defined-inputs to appropriately determine controls
    and properly write fraction infectivity file.
    """
    args = vars(parse_args())
    excelfile = args['infile']
    outfile = args['outfile']
    sheet_map_file = args['sheet_map']
    layouts_dir = args['plate_layouts_dir']

    if not os.path.isfile(excelfile):
        raise ValueError(f"Cannot find `excelfile`{excelfile}.")

    sheet_data = pd.read_excel(excelfile,
                               sheet_name=None,  # read all sheets
                               engine='openpyxl',
                               skiprows=range(0, 22),
                               index_col=0,
                               nrows=8
                               )

    sheet_map_df = pd.read_csv(sheet_map_file)

    required_cols = ['Plate', 'Sample', 'SampleNum', 'Virus', 'PlateLayout',
                     'StartDil', 'DilFactor']
    for col in required_cols:
        if col not in sheet_map_df.columns:
            raise ValueError(f"Required column {col} not in sample map.")

    extra_sheets = set(sheet_data) - set(sheet_map_df['Plate'])
    if extra_sheets:
        raise ValueError(f"`excelfile` {excelfile} has the following extra "
                         f"sheets not in `sheet_mapping`: {extra_sheets}")

    fract_infect_dict = {'serum': [], 'virus': [], 'replicate': [],
                         'concentration': [], 'fraction infectivity': []}

    for plate in sheet_data:
        plate_map = sheet_map_df[sheet_map_df['Plate'] == plate]
        layout_df = pd.read_csv(f"{layouts_dir}/" +
                                f"{plate_map['PlateLayout'].iloc[0]}"
                                ).astype(str)  # in case some cols all number
        plate_df = sheet_data[plate].reset_index(drop=True)
        plate_fract_infect = pd.DataFrame(index=plate_df.index,
                                          columns=plate_df.columns)

        # get background locations
        bg_locs = get_locs(layout_df, 'neg_ctrl')
        if len(bg_locs) == 0:
            if args['allow_cells_only_bg']:
                bg_locs = get_locs(layout_df, 'cells_only')
            else:
                raise ValueError('no "neg_ctrl" background wells')

        bg_rlus = []
        for loc in bg_locs:
            bg_rlus.append(plate_df.at[loc[0], int(loc[1])])

        # get average of bg RLUs and subtract from plate readings
        bg = np.average(bg_rlus)
        plate_bg_sub = plate_df - bg

        # Get locations for positive control (no Ab) wells
        pos_locs = get_locs(layout_df, 'pos_ctrl')
        pos_cols = {loc[1] for loc in pos_locs}
        pos_idxs = {loc[0] for loc in pos_locs}

        # Determine plate orientation from positive control layout
        if len(pos_cols) < len(pos_idxs):
            orientation = 'V'
        elif len(pos_cols) > len(pos_idxs):
            orientation = 'H'
        else:
            raise ValueError("Unable to determine plate orientation from "
                             f"positive control locations ({pos_locs}).")

        # Get values for pos ctrl wells and put in own df
        pos_ctrl_values = pd.DataFrame(index=plate_df.index,
                                       columns=plate_df.columns)
        for pos_loc in pos_locs:
            pos_ctrl_values.at[pos_loc[0], int(pos_loc[1])] = (
                    plate_bg_sub.at[pos_loc[0], int(pos_loc[1])]
                    )

        # Calculate fraction infectivity based on positive ctrl values
        if orientation == 'V':
            pos_ctrl_series = pos_ctrl_values.mean(axis=1, skipna=True)
            for row in plate_fract_infect.index:
                plate_fract_infect.loc[row] = (plate_bg_sub.loc[row] /
                                               pos_ctrl_series[row])

        elif orientation == 'H':
            pos_ctrl_series = pos_ctrl_values.mean(axis=0, skipna=True)
            for col in plate_fract_infect.columns:
                plate_fract_infect[col] = (plate_bg_sub[col] /
                                           pos_ctrl_series[col])

        else:
            raise ValueError(f"Invalid orientation of {orientation}")

        # Get locations for samples and add info to fract_infect_dict
        sample_nums = plate_map['SampleNum'].tolist()
        for sample in sample_nums:
            sample_locs = get_locs(layout_df, str(sample))
            sample_count = len(sample_locs)
            if sample_count <= 0:
                raise ValueError(f"no sample counts for {sample}")
            start_dil = (plate_map[plate_map['SampleNum'] == sample]
                                  ['StartDil'].iloc[0])
            dil_factor = (plate_map[plate_map['SampleNum'] == sample]
                                   ['DilFactor'].iloc[0])

            sample_cols = {loc[1] for loc in sample_locs}
            sample_idxs = {loc[0] for loc in sample_locs}

            # Assuming all of one replicate is either in the same row or col
            reps = min(len(sample_cols), len(sample_idxs))
            if not sample_count % reps == 0:
                raise ValueError("Sample number not evenly divisible by"
                                 f"assumed number of reps {reps}.")

            # Add sample, virus, replicate, and concentration info
            fract_infect_dict['serum'].extend(
                    plate_map[plate_map['SampleNum'] == sample]['Sample']
                    .to_list()*sample_count)
            fract_infect_dict['virus'].extend(
                    plate_map[plate_map['SampleNum'] == sample]['Virus']
                    .to_list()*sample_count)
            for i in range(1, reps+1):
                fract_infect_dict['replicate'].extend(
                       list(itertools.repeat(i, sample_count//reps)))
                fract_infect_dict['concentration'].extend(
                        [(start_dil/(dil_factor**x))
                         for x in range(sample_count//reps)])

            # Add fraction infectivities
            if orientation == 'V':
                for col in sample_cols:
                    fract_infect_dict['fraction infectivity'].extend(
                            plate_fract_infect[int(col)].dropna().to_list())
            elif orientation == 'H':
                for row in sample_idxs:
                    fract_infect_dict['fraction infectivity'].extend(
                            plate_fract_infect.loc[row].dropna().to_list())
            else:
                raise ValueError(f"Invalid orientation: {orientation}")

    assert len(fract_infect_dict['serum']) == \
           len(fract_infect_dict['virus']) == \
           len(fract_infect_dict['replicate']) == \
           len(fract_infect_dict['concentration']) == \
           len(fract_infect_dict['fraction infectivity']), \
           f"Error in making fract_infect_dict:\n{fract_infect_dict}"

    fract_infect_df = pd.DataFrame.from_dict(fract_infect_dict)
    fract_infect_df.to_csv(outfile, float_format='%.4g')


if __name__ == '__main__':
    main()
