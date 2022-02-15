# Uploading FASTQ files to the SRA

Details of how the raw PacBio sequencing files were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
The submitted files are in BioProject [PRJNA770094](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770094).

The Python Jupyter notebook [upload_to_SRA.ipynb](upload_to_SRA.ipynb) has instructions and does the uploading.

Because the FTP upload takes a while, you may want to run the Jupyter notebook using `slurm` so there is no timeout with::

    sbatch --wrap="jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 upload_to_SRA.ipynb" --time 2-0
