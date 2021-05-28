# Deep mutational scanning of the SARS-CoV-2 RBD in variant backgrounds.

Analysis of deep mutational scanning of barcoded codon variants of SARS-CoV-2 RBD in variatn backgrounds.

The 4 backgrounds examined here are:
  1. Wuhan-Hu-1 background ("wildtype")
  2. E484K
  3. N501Y (found in the B.1.1.7 lineage)
  4. K417N-E484K-N501Y (found in the B.1.351 lineage)

Study and analysis by Tyler Starr, Allie Greaney, [Jesse Bloom](https://research.fhcrc.org/bloom/en.html), and co-authors.

## Summary of workflow and results
For a summary of the workflow and links to key results files, [click here](results/summary/summary.md).
Reading this summary is the best way to understand the analysis.

## Running the analysis
The analysis consists of three components, all of which are contained in this repository:

 1. Instructions to build the computing environment.

 2. The computer code itself.

 3. The required input data.

First, set up the computing environment, which is partially done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
The environment is specified in [environment.yml](environment.yml).
If you have not previously built the conda environment, then build the environment specified in [environment.yml](environment.yml) to `./env` with:

    conda env create -f environment.yml -p ./env

Then activate it with:

    conda activate ./env

If you've previously built the environment into `./env`, just do the activation step.

Setting up the `conda` environment above installs everything to run all parts of the analysis **except** the `R` markdown notebooks.
For those, the pipeline currently uses the Fred Hutch computing cluster module `R/3.6.2-foss-2019b` as specified in `Snakefile`.
That module is not packaged with this repo, so if you aren't on the Fred Hutch cluster you'll have to create a similar `R` environment yourself (all the `R` packages are listed at the beginning of their output in the [summary results](results/summary/summary.md).

Now you can run the entire analysis.
The analysis consists primarily of a series of Jupyter notebooks and R markdown in to the top-level directory along with some additional code in [Snakefile](Snakefile).
You can run the analysis by using [Snakemake](https://snakemake.readthedocs.io) to run [Snakefile](Snakefile), specifying the conda environment in `./env`, as in:

    snakemake --use-conda --conda-prefix ./env

However, you probably want to use the server to help with computationally intensive parts of the analysis.
To run using the cluster configuration for the Fred Hutch server, simply run the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash), which executes [Snakefile](Snakefile) in a way that takes advantage of the Hutch server resources.
This bash script also automates the environment building steps above, so really all you have to do is run this script.
You likely want to submit [run_Hutch_cluster.bash](run_Hutch_cluster.bash) itself to the cluster (since it takes a while to run) with:

    sbatch -t 7-0 run_Hutch_cluster.bash

### Configure `.git` to not track Jupyter notebook metadata
To simplify git tracking of Jupyter notebooks, we have added the filter described [here](https://stackoverflow.com/questions/28908319/how-to-clear-an-ipython-notebooks-output-in-all-cells-from-the-linux-terminal/58004619#58004619) to strip notebook metadata to [.gitattributes](.gitattributes) and [.gitconfig](.gitconfig).
The **first time** you check out this repo, run the following command to use this configuration (see [here](https://stackoverflow.com/a/18330114)):

   git config --local include.path ../.gitconfig

Then don't worry about it anymore.

## Configuring the analysis
The configuration for the analysis is specifed in [config.yaml](config.yaml).
This file defines key variables for the analysis, and should be relatively self-explanatory.
You should modify the analysis by changing this configuration file; do **not** hard-code crucial experiment-specific variables within the Jupyter notebooks or `Snakefile`.

The input files pointed to by [config.yaml](config.yaml) are in the [./data/](data) subdirectory.
See the [./data/README.md](./data/README.md) file for details.

Note that the raw sequencing data are on the SRA in BioProject [PRJNA639956](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA639956) as well as on the Hutch cluster.

## Cluster configuration
There is a cluster configuration file [cluster.yaml](cluster.yaml) that configures [Snakefile](Snakefile) for the Fred Hutch cluster, as recommended by the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
The [run_Hutch_cluster.bash](run_Hutch_cluster.bash) script uses this configuration to run [Snakefile](Snakefile).
If you are using a different cluster than the Fred Hutch one, you may need to modify the cluster configuration file.

## Notebooks that perform the analysis
The Jupyter notebooks and R markdown dscripts that perform most of the analysis are in this top-level directory with the extension `*.ipynb` or `*.Rmd`.
These notebooks read the key configuration values from [config.yaml](config.yaml).

There is also a [./scripts/](scripts) subdirectory with related scripts.

The notebooks need to be run in the order described in [the workflow and results summary](results/summary/summary.md).
This will occur automatically if you run them via [Snakefile](Snakefile) as described above.

## Results
Results are placed in the [./results/](results) subdirectory.
Many of the files created in this subdirectory are not tracked in the `git` repo as they are very large.
However, key results files are tracked as well as a summary that shows the code and results.
Click [here](./results/summary/summary.md) to see that summary.

The large results files are tracked via [git-lfs](https://git-lfs.github.com/).
This requires `git-lfs` to be installed, which it is in the `conda` environment specified by [environment.yml](environment.yml).
The following commands were then run:

    git lfs install

You may need to run this if you are tracking these files and haven't installed `git-lfs` in your user account.
Then the large results files were added for tracking with: xxx

## Updating the conda environment
The [environment.yml](environment.yml) file contains a fully pinned conda environment.
An environment without all of the versions pinned is in [environment_unpinned.yml](environment_unpinned.yml).
If you need to update the environment, the suggested way to do it is add the new requirement to [environment_unpinned.yml](environment_unpinned.yml), then build that environment and finally export the pinned version:

    conda env create -f environment_unpinned.yml --prefix ./env

Then activate the environment with:

    conda activate ./env

Finally, export the pinned version with:

    conda env export --prefix ./env > environment.yml
