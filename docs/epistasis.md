---
layout: epistasis
permalink: /epistatic-shifts/
---

## Overview 

You can use this tool to explore epistatic shifts in mutational effects on ACE2-binding affinity (-log10 Kd) between SARS-CoV-2 receptor-binding domain (RBD) variants. Once you've selected a comparison of interest, you can investigate the mutation-level epistatic shifts. 

### Instructions

To use this tool, select two SARS-CoV-2 variants that you wish to compare between. To do this, simply select a *'comparator'* in the dropdown menu below the plot and select a *'variant'* by clicking on a variant name in the legend above the plot. Now you're comparing the epistatic shift (Jensen-Shannon divergence between ACE2 binding affinities) at each RBD site between the variant backgrounds you selected. 

Now that you've selected two variants to compare between, you can investigate site level differences by clicking on the points in the higlighted line plot. Simply click on a point – you should see the size of the point change indicating your selection – then you will see the differences in affinities of each of the 20 amino acids measured in each variant appear in the scatter plot on the right. Hover over individual mutations to see exact numerical details.

### Technical Details

The epistatic shift is calculated as the Jensen-Shannon divergence in the set of Boltzmann-weighted affinities for all amino acids at each site. Mutation affinities were experimentally measured via high-throughput ACE2-binding titrations with yeast-displayed RBDs.

More information can be found in our preprint [here](LINK NOT YET AVAILABLE). Raw data can be found [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants/blob/main/results/epistatic_shifts/JSD_versus_Wuhan1_by_target.csv) for a table of all pairwise RBD epistatic shifts, and [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants/blob/main/results/final_variant_scores/final_variant_scores.csv) for individual measurements of RBD mutant affinities.
