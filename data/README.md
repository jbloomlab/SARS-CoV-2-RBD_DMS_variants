# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - `PacBio_amplicons_*_.gb`: the amplicons being sequenced by PacBio.
     There are 4 different backgrounds: Wuhan-1, E484K, N501Y, and B.1.351 (which has a K417N-E484K-N501Y triple mutation)

   - `feature_parse_specs_*_.yaml`: how to parse the amplicon when handling the PacBio data.

   - [PacBio_runs.csv](PacBio_runs.csv): list of the PacBio runs used to call the variants.

   - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples.

## Structures

   - [6m0j.pdb](6m0j.pdb): ACE2-bound SARS-CoV-2 RBD structure
   
   - [./VOC_structures](./VOC_structures): structurally aligned RBDs of ACE2-bound RBD structures for the Wuhan-Hu-1, alpha, beta, and omicron structures. See the description in the [structural_shifts](../structural_shifts.Rmd) analysis of how these files were created.

   