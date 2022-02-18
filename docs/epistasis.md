---
layout: epistasis
permalink: /epistatic-shifts/
---

## Overview 

You can use this tool to explore the epistatic shift in ACE2 binding affinity (-log10 Kd) between the RBD libraries of your choosing. Once you've selected a comparison of interest, you can investigate the site level differences. 

### Instructions

To use this tool, select two SARS-CoV-2 variants that you wish to compare between. To do this, simply select a *'comparator'* in the dropdown menu below the plot and select a *'variant'* by clicking on a variant name in the legend above the plot. Now you're comparing the epistatic shift (Jensen-Shannon Distance between ACE2 binding affinities) between the variant backgrounds you selected. 

Now that you've selected two variants to compare between, you can investigate site level differences by clicking on the points in the higlighted line plot. Simply click on a point – you should see the size of the point change indicating your selection – then you will see the amino acid level differences at that site appear in the scatter plot on the left. 

### Technical Details

The epistatic shift is calculated as the Jensen-Shannon divergence in the set of Boltzmann-weighted affinities for all amino acids at each site.