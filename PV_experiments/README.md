# Spike-pseudotyped lentiviral entry and neutralization experimental data

Analysis of validation data outside the primary computational pipeline. 

Rescued eight viruses (plus no VEP controls) in triplicate. Used one set (N501Y/Y449H cycle) for both titering on ACE2-high and ACE2-low expressing cell lines, and neutralization by a panel of monoclonal antibodies (three expected to be escaped by Y449H to various degrees from prior deep mutational scanning experiments, and one negative control [S2E12] that is not expacted to be escaped). Used the other set (R498Q/Y501N reversion cycle in omicron) for titering experiments.

Raw data on entry titers is given in the [titer_data](./titer_data) subdirectory. This data is analyzed in the [pseudovirus_titer.Rmd](./pseudovirus_titer.Rmd) script, with results in the [results](./results) directory.

Raw data on antibody neutralization is given in the [neut_data](./neut_data) subdirectory, and neutralization curves are fit by the [analyze_neut_data.ipynb](./analyze_neut_data.ipynb) notebook, with results in the [results/neut_titers](./results/neut_titers) subdirectory.