# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import glob
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.2.4
    Using dms_variants version 0.8.9


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_Wuhan_Hu_1'], na_filter=None)
variants = variants.append(pd.read_csv(config['codon_variant_table_file_E484K'], na_filter=None))
variants = variants.append(pd.read_csv(config['codon_variant_table_file_N501Y'], na_filter=None))
variants = variants.append(pd.read_csv(config['codon_variant_table_file_B1351'], na_filter=None))

variants = variants.reset_index(drop=True)

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAAAATTTAA</td>
      <td>4</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAAACGCGTA</td>
      <td>3</td>
      <td>GAA154ACT</td>
      <td>E154T</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAAACTCCAA</td>
      <td>2</td>
      <td>TTT156ATG</td>
      <td>F156M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAACCGATTA</td>
      <td>2</td>
      <td>CAG84GAA</td>
      <td>Q84E</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAACGGATGA</td>
      <td>1</td>
      <td>GCT14GGT</td>
      <td>A14G</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAGAGACAATTCGTT</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAAGCCGGATTCGTAC</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAATATGAAAGATACA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACAGCCGATTTACAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAGAGCATAAGCCCCA</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 347 duplicated barcodes.Started with 330417 barcodes:
    After removing duplicates, there are 329723 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_Wuhan_Hu_1'],
                                     feature_parse_specs=config['feature_parse_specs_Wuhan_Hu_1'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files expanded by glob
barcode_runs = (
    pd.read_csv(config['barcode_runs'])
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            )
    )

display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1_prefix</th>
      <th>R1_samplename</th>
      <th>R1_postfix</th>
      <th>R1</th>
      <th>n_R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>210816</td>
      <td>1410589</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s01-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s01-b1_S1_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>210816</td>
      <td>533374</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s01-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s01-b2_S2_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>210816</td>
      <td>1697538</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s01-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s01-b3_S3_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>210816</td>
      <td>6497928</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s01-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s01-b4_S4_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>210816</td>
      <td>1587111</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s02-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s02-b1_S5_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>210816</td>
      <td>882688</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s02-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s02-b2_S6_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>210816</td>
      <td>2168488</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s02-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s02-b3_S7_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>210816</td>
      <td>5908495</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s02-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s02-b4_S8_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>210816</td>
      <td>2236680</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s03-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s03-b1_S9_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>210816</td>
      <td>1125722</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s03-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s03-b2_S10_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>210816</td>
      <td>3404113</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s03-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s03-b3_S11_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>210816</td>
      <td>3354021</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s03-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s03-b4_S12_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>210816</td>
      <td>3228357</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s04-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s04-b1_S13_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>210816</td>
      <td>3777633</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s04-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s04-b2_S14_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>210816</td>
      <td>2507212</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s04-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s04-b3_S15_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>210816</td>
      <td>522258</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s04-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s04-b4_S16_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>210816</td>
      <td>5380520</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s05-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s05-b1_S17_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>210816</td>
      <td>4516641</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s05-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s05-b2_S18_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>210816</td>
      <td>458022</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s05-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s05-b3_S19_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>210816</td>
      <td>2412</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s05-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s05-b4_S20_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>210816</td>
      <td>6953452</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s06-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s06-b1_S21_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>210816</td>
      <td>3350686</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s06-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s06-b2_S22_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>210816</td>
      <td>4212</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s06-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s06-b3_S23_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>210816</td>
      <td>304</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s06-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s06-b4_S24_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>210816</td>
      <td>8037439</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s07-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s07-b1_S25_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>210816</td>
      <td>2440419</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s07-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s07-b2_S26_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>210816</td>
      <td>1587</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s07-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s07-b3_S27_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>210816</td>
      <td>470</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s07-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s07-b4_S28_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>210816</td>
      <td>8013743</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s08-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s08-b1_S29_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>210816</td>
      <td>2313562</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s08-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s08-b2_S30_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>210816</td>
      <td>1458</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s08-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s08-b3_S31_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>210816</td>
      <td>276</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s08-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s08-b4_S32_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>210816</td>
      <td>8191465</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s09-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s09-b1_S33_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>210816</td>
      <td>2003725</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s09-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s09-b2_S34_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>210816</td>
      <td>882</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s09-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s09-b3_S35_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>210816</td>
      <td>275</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s09-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s09-b4_S36_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>210816</td>
      <td>747246</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s10-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s10-b1_S37_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>210816</td>
      <td>632407</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s10-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s10-b2_S38_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>210816</td>
      <td>1569589</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s10-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s10-b3_S39_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>210816</td>
      <td>7076796</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s10-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s10-b4_S40_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>210816</td>
      <td>508269</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s11-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s11-b1_S41_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>210816</td>
      <td>1360633</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s11-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s11-b2_S42_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>210816</td>
      <td>2122153</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s11-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s11-b3_S43_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>210816</td>
      <td>6133995</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s11-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s11-b4_S44_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>210816</td>
      <td>1471208</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s12-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s12-b1_S45_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>210816</td>
      <td>1351829</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s12-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s12-b2_S46_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>210816</td>
      <td>3465388</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s12-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s12-b3_S47_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>210816</td>
      <td>3785496</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s12-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s12-b4_S48_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>210816</td>
      <td>2711471</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s13-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s13-b1_S49_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>210816</td>
      <td>3855154</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s13-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s13-b2_S50_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>210816</td>
      <td>3017813</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s13-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s13-b3_S51_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>210816</td>
      <td>564315</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s13-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s13-b4_S52_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>210816</td>
      <td>6181987</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s14-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s14-b1_S53_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>210816</td>
      <td>3991649</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s14-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s14-b2_S54_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>210816</td>
      <td>438596</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s14-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s14-b3_S55_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>210816</td>
      <td>1844</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s14-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s14-b4_S56_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>210816</td>
      <td>8575711</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s15-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s15-b1_S57_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>210816</td>
      <td>1572345</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s15-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s15-b2_S58_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>210816</td>
      <td>2862</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s15-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s15-b3_S59_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>210816</td>
      <td>190</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s15-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s15-b4_S60_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>210816</td>
      <td>9233079</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s16-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s16-b1_S61_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>210816</td>
      <td>870320</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s16-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s16-b2_S62_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>210816</td>
      <td>617</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s16-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s16-b3_S63_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>210816</td>
      <td>129</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s16-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s16-b4_S64_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>210816</td>
      <td>9320303</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s17-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s17-b1_S65_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>210816</td>
      <td>772536</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s17-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s17-b2_S66_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>210816</td>
      <td>571</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s17-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s17-b3_S67_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>210816</td>
      <td>141</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s17-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s17-b4_S68_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>210816</td>
      <td>9320303</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s18-b1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s18-b1_S69_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>210816</td>
      <td>772536</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s18-b2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s18-b2_S70_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>210816</td>
      <td>571</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s18-b3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s18-b3_S71_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>210816</td>
      <td>141</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210816_s18-b4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210816_s18-b4_S72_R1_001.fastq.gz]</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>210811</td>
      <td>1480000</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib1_bin1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin1_2_S74_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin1_1_S73_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin1_3_S75_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>210811</td>
      <td>3840000</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib1_bin2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin2_1_S76_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin2_2_S77_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin2_3_S78_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>210811</td>
      <td>3272500</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib1_bin3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin3_2_S80_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin3_1_S79_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin3_3_S81_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>210811</td>
      <td>3336000</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib1_bin4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin4_3_S84_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin4_1_S82_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib1_bin4_2_S83_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>210811</td>
      <td>1280000</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib2_bin1</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin1_2_S86_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin1_1_S85_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin1_3_S87_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>210811</td>
      <td>3120000</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib2_bin2</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin2_3_S90_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin2_1_S88_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin2_2_S89_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>210811</td>
      <td>3025000</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib2_bin3</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin3_2_S92_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin3_1_S91_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin3_3_S93_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>210811</td>
      <td>3360000</td>
      <td>/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/</td>
      <td>210811_lib2_bin4</td>
      <td>*R1*.fastq.gz</td>
      <td>[/shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin4_2_S95_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin4_1_S94_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/210903_D00300_1313_AHMJJTBCX3/Unaligned/Project_tstarr/210811_lib2_bin4_3_S96_R1_001.fastq.gz]</td>
      <td>3</td>
    </tr>
  </tbody>
</table>


Make sure all samples have at least one R1 file


```python
if any(barcode_runs['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs.query('n_R1 < 1')}")
```

Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>174804</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>154919</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 8 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>TCAGCACTACGGAAAC</td>
      <td>504</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>ATACTTATGTATAGAC</td>
      <td>498</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CGTCCTTGCTGTCGAG</td>
      <td>480</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>AAGGAAGTAGCCCCTT</td>
      <td>477</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CTACTCTATATTCAAA</td>
      <td>462</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>1401160</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>313727</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>66748</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>27326</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>26207</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="40" valign="top">pool1</th>
      <th>SortSeq_bin1</th>
      <td>25110</td>
      <td>435884</td>
      <td>59972</td>
      <td>33793</td>
      <td>1078183</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>57437</td>
      <td>657901</td>
      <td>142246</td>
      <td>68207</td>
      <td>2936893</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>55814</td>
      <td>648207</td>
      <td>138007</td>
      <td>59785</td>
      <td>2867027</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>59361</td>
      <td>686860</td>
      <td>147224</td>
      <td>65302</td>
      <td>3009319</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>27326</td>
      <td>313727</td>
      <td>66748</td>
      <td>26207</td>
      <td>1401160</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>8411</td>
      <td>92074</td>
      <td>19583</td>
      <td>10902</td>
      <td>422482</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>29568</td>
      <td>338575</td>
      <td>73228</td>
      <td>26256</td>
      <td>1535028</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>114278</td>
      <td>1283460</td>
      <td>287019</td>
      <td>97643</td>
      <td>5846987</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>26390</td>
      <td>303590</td>
      <td>64761</td>
      <td>27399</td>
      <td>1344544</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>20086</td>
      <td>208274</td>
      <td>48355</td>
      <td>19077</td>
      <td>954321</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>32582</td>
      <td>377600</td>
      <td>81080</td>
      <td>35233</td>
      <td>1701203</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>101100</td>
      <td>1130551</td>
      <td>245215</td>
      <td>95102</td>
      <td>5220928</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>37474</td>
      <td>426421</td>
      <td>95264</td>
      <td>37008</td>
      <td>1918414</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>20872</td>
      <td>229121</td>
      <td>48470</td>
      <td>19744</td>
      <td>1032300</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>58421</td>
      <td>681964</td>
      <td>144471</td>
      <td>53352</td>
      <td>2892343</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>58170</td>
      <td>619463</td>
      <td>143851</td>
      <td>58407</td>
      <td>3015276</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>54712</td>
      <td>613188</td>
      <td>125983</td>
      <td>48339</td>
      <td>2643005</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>61637</td>
      <td>741582</td>
      <td>149624</td>
      <td>52777</td>
      <td>3078492</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>40438</td>
      <td>416432</td>
      <td>98779</td>
      <td>34315</td>
      <td>2099291</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>17196</td>
      <td>160562</td>
      <td>43466</td>
      <td>15760</td>
      <td>904630</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>86583</td>
      <td>1029280</td>
      <td>210981</td>
      <td>72204</td>
      <td>4438344</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>75598</td>
      <td>815752</td>
      <td>179072</td>
      <td>69650</td>
      <td>3809744</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>8321</td>
      <td>77775</td>
      <td>19291</td>
      <td>10030</td>
      <td>434653</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>520</td>
      <td>2165</td>
      <td>432</td>
      <td>5786</td>
      <td>8481</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>132843</td>
      <td>1506848</td>
      <td>326909</td>
      <td>111102</td>
      <td>6751350</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>60642</td>
      <td>643630</td>
      <td>143512</td>
      <td>50229</td>
      <td>3038763</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>325</td>
      <td>2291</td>
      <td>549</td>
      <td>619</td>
      <td>11284</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>395</td>
      <td>1774</td>
      <td>374</td>
      <td>3767</td>
      <td>7071</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>152043</td>
      <td>1664791</td>
      <td>368249</td>
      <td>127097</td>
      <td>7612770</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>42543</td>
      <td>469731</td>
      <td>102926</td>
      <td>42405</td>
      <td>2140046</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>212</td>
      <td>504</td>
      <td>112</td>
      <td>527</td>
      <td>1702</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>168</td>
      <td>926</td>
      <td>218</td>
      <td>407</td>
      <td>4318</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>133313</td>
      <td>1561482</td>
      <td>335070</td>
      <td>113355</td>
      <td>7167389</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>45834</td>
      <td>495828</td>
      <td>106511</td>
      <td>46247</td>
      <td>2234706</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>110</td>
      <td>475</td>
      <td>76</td>
      <td>999</td>
      <td>1395</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>189</td>
      <td>1342</td>
      <td>73</td>
      <td>2133</td>
      <td>282</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>156059</td>
      <td>1734792</td>
      <td>374595</td>
      <td>129919</td>
      <td>7906099</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>31746</td>
      <td>358803</td>
      <td>78777</td>
      <td>26975</td>
      <td>1619554</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>323</td>
      <td>472</td>
      <td>111</td>
      <td>2812</td>
      <td>1796</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>58</td>
      <td>223</td>
      <td>29</td>
      <td>779</td>
      <td>495</td>
    </tr>
    <tr>
      <th rowspan="40" valign="top">pool2</th>
      <th>SortSeq_bin1</th>
      <td>21357</td>
      <td>416194</td>
      <td>50353</td>
      <td>56206</td>
      <td>841152</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>55516</td>
      <td>722369</td>
      <td>135620</td>
      <td>64611</td>
      <td>2670367</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>52354</td>
      <td>636939</td>
      <td>129546</td>
      <td>54545</td>
      <td>2498610</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>60506</td>
      <td>733364</td>
      <td>147327</td>
      <td>68770</td>
      <td>2950412</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>9148</td>
      <td>116061</td>
      <td>22174</td>
      <td>8739</td>
      <td>440940</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>8361</td>
      <td>109375</td>
      <td>21010</td>
      <td>7065</td>
      <td>417275</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>26059</td>
      <td>324493</td>
      <td>62823</td>
      <td>21568</td>
      <td>1269043</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>118677</td>
      <td>1489368</td>
      <td>289294</td>
      <td>100511</td>
      <td>5935478</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>6748</td>
      <td>89409</td>
      <td>17275</td>
      <td>6191</td>
      <td>346367</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>24962</td>
      <td>317080</td>
      <td>59832</td>
      <td>20368</td>
      <td>1203344</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>35127</td>
      <td>447118</td>
      <td>84464</td>
      <td>29032</td>
      <td>1756846</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>107767</td>
      <td>1323046</td>
      <td>265812</td>
      <td>87325</td>
      <td>5231867</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>29241</td>
      <td>377070</td>
      <td>72161</td>
      <td>24281</td>
      <td>1440718</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>19015</td>
      <td>257848</td>
      <td>47223</td>
      <td>16718</td>
      <td>984381</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>55913</td>
      <td>681297</td>
      <td>133706</td>
      <td>48191</td>
      <td>2780367</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>59244</td>
      <td>750270</td>
      <td>142794</td>
      <td>47724</td>
      <td>2902784</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>39070</td>
      <td>493816</td>
      <td>93395</td>
      <td>33800</td>
      <td>1917323</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>65058</td>
      <td>781755</td>
      <td>156347</td>
      <td>53691</td>
      <td>3204708</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>44244</td>
      <td>543624</td>
      <td>106731</td>
      <td>46589</td>
      <td>2159414</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>10297</td>
      <td>140523</td>
      <td>24411</td>
      <td>10551</td>
      <td>489256</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>64870</td>
      <td>811995</td>
      <td>162398</td>
      <td>52666</td>
      <td>3263911</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>68241</td>
      <td>846107</td>
      <td>166988</td>
      <td>61202</td>
      <td>3335463</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>6694</td>
      <td>99010</td>
      <td>17296</td>
      <td>7772</td>
      <td>335456</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>271</td>
      <td>1398</td>
      <td>217</td>
      <td>3319</td>
      <td>4053</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>144029</td>
      <td>1820090</td>
      <td>360142</td>
      <td>116931</td>
      <td>7281395</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>24716</td>
      <td>328070</td>
      <td>61639</td>
      <td>22814</td>
      <td>1225081</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>135</td>
      <td>294</td>
      <td>41</td>
      <td>3841</td>
      <td>649</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>192</td>
      <td>38</td>
      <td>25</td>
      <td>1824</td>
      <td>124</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>165255</td>
      <td>2099705</td>
      <td>397628</td>
      <td>132052</td>
      <td>8226207</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>16478</td>
      <td>208095</td>
      <td>39972</td>
      <td>16660</td>
      <td>811923</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>119</td>
      <td>518</td>
      <td>107</td>
      <td>1832</td>
      <td>2060</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>1368</td>
      <td>140</td>
      <td>14</td>
      <td>6604</td>
      <td>133</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>142992</td>
      <td>1782858</td>
      <td>347636</td>
      <td>127355</td>
      <td>7073711</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>12721</td>
      <td>165049</td>
      <td>31934</td>
      <td>11120</td>
      <td>635735</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>100</td>
      <td>227</td>
      <td>54</td>
      <td>2175</td>
      <td>936</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>102</td>
      <td>49</td>
      <td>9</td>
      <td>1140</td>
      <td>45</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>149118</td>
      <td>1868460</td>
      <td>371898</td>
      <td>129325</td>
      <td>7373815</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>12078</td>
      <td>156337</td>
      <td>29807</td>
      <td>12012</td>
      <td>605846</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>449</td>
      <td>3042</td>
      <td>227</td>
      <td>5736</td>
      <td>2754</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>133</td>
      <td>121</td>
      <td>18</td>
      <td>2835</td>
      <td>88</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_44_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
