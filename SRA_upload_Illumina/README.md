# Uploading FASTQ files to the SRA as a new BioSample on an existing BioProject

Before publishing your study, you need to make the raw data available.
Here is how to upload FASTQ files to the NCBI [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
The files are uploaded to a recently created BioProject [PRJNA770094](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770094), which is set up to store all deep mutational scanning and mutational antigenic profiling experiments that use the variant RBD SSM libraries described in this work.

These instructions are for uploading Illumina barcode sequencing for the expression Sort-seq and ACE2-binding Tite-seq experiments.

For each new deep mutational scanning or mutational antigenic profiling study with these libraries, you can upload the FASTQ files as a new *BioSample*.
Here are the steps:

## Create the BioSample
Go to the [SRA Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/) and login.

Then under the `Start a new submission` banner, click *BioSample* to create a new *BioSample*.
You are now on a page that has a blue `New Submission` button, which you should select.
This will bring you to a page with some Submitter information; check that it is correct and then hit `Continue`.
You will then be on the `General Information` page.
Select when to release the submission to the public (generally, immediately is OK).
Then click that you are uploading a *Single BioSample* and hit `Continue`.
You will now be at the *Sample Type* page, and you have to select the package that best describes the submission.
Click *Microbe*--although we are studying a pathogen it's not a direct clinical sample but an experiment using yeast on a pathogen, which is why we choose this sample type.
Then click `Continue`.
Now you enter the sample attributes.
For the sample name, provide a short name that describes the sample, such as `variant_RBD_DMS`.
Also provide the rest of the information:

  - Organism: Severe acute respiratory syndrome-related coronavirus

  - strain: various

  - isolation source: plasmid

  - collection date: 2021

  - geographic location: USA

  - sample type: plasmid

Then hit `Continue`.
You will now be on the page to specify the BioProject.
We are adding to an existing BioProject, so enter [PRJNA770094](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770094) as the *Existing BioProject* and hit `Continue`.
Finally, add a sample title, such as "Illumina barcode sequencing for SARS2 variant RBD ACE2 and expression DMS."
Then hit `Continue`, make sure everything looks correct, then hit `Submit`.

After a brief bit of processing, the *BioSample* submission should show up, along with a sample accession that will be in the format of *SAMN25944367*.
Add this sample accession to [upload_config.yaml](upload_config.yaml) as the value for the *biosample_accession* key.

## Upload the sequencing data
Again go to the [SRA Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/) and login.
This time, under the `Start a new submission` banner, click *Sequence Read Archive* to upload the actual sequencing data.

You are now on a page with a `New submission` button, which you should click.
Check that the submitter information is correct, then click `Continue`.
You will now be at the *General Information* page.
We are adding to an existing BioProject, so enter [PRJNA770094](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770094) as the *Existing BioProject* and hit `Continue`.
Then for the question of whether you already registered a *BioSample*, also select yes.
Then select a *Release Date* depending on whether you want to release the results immediately (usually fine) or at some future date.
Then click `Continue`.

You will now be at a page that asks how you want to provide the SRA Metadata.
Click to *Upload a file*.
The next section describes how to create this table.

## Create SRA metadata submission table
In order to create the submission table, you need to specify which Illumina barcode runs we are uploading and what the sample name should be.
This can be done by editing the [upload_config.yaml](upload_config.yaml) file.
This YAML has two entries:

 - The `barcode_runs` file, which is a CSV in the format of the repo's master barcode run file at [../data/barcode_runs.csv](../data/barcode_runs.csv). But importantly, this file should **only** have the barcode runs for the samples of interest for this *BioSample*. If you have already made a subset repo (see [../subset_data](../subset_data)) that just has those samples of interest, then you can put the path to the master [../data/barcode_runs.csv](../data/barcode_runs.csv) file for the whole repo. Otherwise you want to create a CSV file that just holds the barcode runs of interest (this can be done as described in [../subset_data](subset_data)) and put the path to that.

 - The `sample_name` which should be a concise description (no spaces) for the *BioSample*. This should be the same name used when creating the *BioSample* above.

After configuring [upload_config.yaml](upload_config.yaml), simply run the Jupyter notebook [create_submission_sheet.ipynb](create_submission_sheet.ipynb), which will create two output files:

  - [SRA_submission_spreedsheet.tsv](SRA_submission_spreedsheet.tsv) which is for uploading to the SRA.

  - [FASTQs_to_upload.csv](FASTQs_to_upload.csv), which lists all the FASTQs that need to be uploaded.

## Upload the submission table
Now return to the [SRA Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/sra/) webpage, and if needed navigate back to your submission.
You should still be at the `BioSample Attributes` step, and see a `Choose file` box to upload your BioSample attributes.
Use this click box to upload the [SRA_submission_spreedsheet.tsv](SRA_submission_spreedsheet.tsv) file that you just created in the step above.
Then click `Continue`.

After a little while, you should now get a page that asks you how you want to upload the files for this submission.
Click the option for *FTP or Aspera Command Line file preload*.

If you click on the `+` FTP upload instructions, you will see details.
Add the `Username` and `account folder` provided in these instructions to [upload_config.yaml](upload_config.yaml) as the values for the *ftp_username* and *ftp_account_folder* keys.
Also add a value for the *ftp_subfolder* that is meaningful for this particular submission, such as *variant_RBD_Illumina_bcs*.
Finally, put the FTP password as plain text in a file called `ftp_password.txt` which is **not** tracked in this repo for privacy.

## Upload the sequencing data
Now we need to upload the actual sequencing data.
This is done by the Jupyter notebook [make_and_upload_tar.ipynb](make_and_upload_tar.ipynb).
This notebook should run completely without error if you have set the [upload_config.yaml](upload_config.yaml) and `ftp_password.txt` files properly as described above.
It first creates a very large `*.tar` file called `SRA_submission.tar` that contains all the FASTQs.
It then uses FTP to upload them to the SRA.
You can run the notebook interactively, but it will take a little while so make sure it doesn't time out.
If you want to instead submit it via `slurm`, do this with:

    sbatch --wrap="jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 make_and_upload_tar.ipynb" --time 1-0

After you have finished running [make_and_upload_tar.ipynb](make_and_upload_tar.ipynb), check carefully to make sure the FTP upload was completed.
If needed, you can manually log into the FTP site to see the file and use `ls` to see the size of what has been transferred.

Finally, return to the SRA submission webpage for the reads, and check the blue `Select preload folder` box.
Note that you need to wait about 10 minutes for the pre-load folder to become visible.
The click to select the folder you created (this is the *ftp_subfolder* defined in [upload_config.yaml](upload_config.yaml)) and click *Use selected folder*.
Finally, check `Autofinish submission` box and hit `Continue`.
You will get a warning that files are missing since you uploaded a `*.tar` archive; do **not** worry about this and just click `Continue` again.
The webpage will then indicate it is extracting files from the `*.tar`, so wait for this to finish.
It should then show that your submission is complete and just waiting for processing.

You then probably want to delete the `SRA_submission.tar` file as it is very large.

