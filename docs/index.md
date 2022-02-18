---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
---

## Overview 

SARS-CoV-2 has evolved variants with substitutions in the spike receptor-binding domain (RBD) that impact its affinity for ACE2 receptor and recognition by antibodies. These substitutions could also shape future evolution by modulating the effects of mutations at other sites—a phenomenon called epistasis. To investigate this possibility, we performed deep mutational scans to measure the impact on ACE2 binding of all single amino-acid mutations in the Wuhan-Hu-1, Alpha, Beta, Delta, and Eta variant RBDs.

Here, we link to two interactive visualizations that help explore the results of these deep mutation scans and the epistatic shifts between variants. The manuscript detailing the results of these experiments is published [here](). 

### Instructions 

We have made two tools to help visualize the data from our deep mutations scans:

1. A set of interactive heatmaps that you can use to explore the change in binding affinity ($$\Delta$$log10 $$K_D$$) to ACE2 or the change in expression of the RBD (log10(MFI)) between the mutant and wildtype amino acids from each experiment. To use this tool, click [here]({{ site.baseurl }}{% link heatmaps.md %}).

2. An interactive widget that you can use to visualize the epistatic shift in ACE2 binding affinity (-log10 $$K_D$$) between the five variant RBD backgrounds. To use this tool, click [here]({{ site.baseurl }}{% link epistasis.md %}).  

### Data

If you are interested in the raw data from our study, you can find the ACE2 binding affinity (-log10 $$K_D$$) and RBD expression (log10(MFI)) for each experiment [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants/blob/visualizations/results/final_variant_scores/final_variant_scores.csv). You can find the data used to plot the epistatic shifts between variant backgrounds [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants/blob/visualizations/results/epistatic_shifts/JSD_versus_Wuhan1_by_target.csv). 
