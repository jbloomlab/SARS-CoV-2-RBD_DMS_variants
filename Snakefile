"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['bg'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# barcode runs with R1 files expanded by glob
barcode_runs_expandR1 = (
    barcode_runs
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

if any(barcode_runs_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_expandR1.query('n_R1 < 1')}")    

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        get_mut_bind_expr=config['mut_bind_expr'],
        process_ccs_Wuhan_Hu_1=nb_markdown('process_ccs_Wuhan_Hu_1.ipynb'),
        process_ccs_E484K=nb_markdown('process_ccs_E484K.ipynb'),
        process_ccs_N501Y=nb_markdown('process_ccs_N501Y.ipynb'),
        process_ccs_B1351=nb_markdown('process_ccs_B1351.ipynb'),
        barcode_variant_table_Wuhan_Hu_1=config['codon_variant_table_file_Wuhan_Hu_1'],
        barcode_variant_table_E484K=config['codon_variant_table_file_E484K'],
        barcode_variant_table_N501Y=config['codon_variant_table_file_N501Y'],
        barcode_variant_table_B1351=config['codon_variant_table_file_B1351'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        fit_titrations='results/summary/compute_binding_Kd.md',
        variant_Kds_file=config['Titeseq_Kds_file'],
        calculate_expression='results/summary/compute_expression_meanF.md',
        variant_expression_file=config['expression_sortseq_file'],
        collapse_scores='results/summary/collapse_scores.md',
        mut_phenos_file=config['final_variant_scores_mut_file'],
        UShER_tree=config['UShER_tree'],
        refseq=config['UShER_ref']
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:

            1. Get prior Wuhan-1 RBD DMS mutation-level [binding and expression data]({path(input.get_mut_bind_expr)}). 
            
            2. Process PacBio CCSs for each background: [Wuhan_Hu_1]({path(input.process_ccs_Wuhan_Hu_1)}), [E484K]({path(input.process_ccs_E484K)}), [N501Y]({path(input.process_ccs_N501Y)}), [B.1.351]({path(input.process_ccs_B1351)}). Creates barcode-variant lookup tables for each background: [Wuhan_Hu_1]({path(input.barcode_variant_table_Wuhan_Hu_1)}), [E484K]({path(input.barcode_variant_table_E484K)}), [N501Y]({path(input.barcode_variant_table_N501Y)}), [B.1.351]({path(input.barcode_variant_table_B1351)}).
            
            3. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            4. [Fit titration curves]({path(input.fit_titrations)}) to calculate per-barcode K<sub>D</sub>, recorded in [this file]({path(input.variant_Kds_file)}).
            
            5. [Analyze Sort-seq]({path(input.calculate_expression)}) to calculate per-barcode RBD expression, recorded in [this file]({path(input.variant_expression_file)}).
            
            6. [Derive final genotype-level phenotypes from replicate barcoded sequences]({path(input.collapse_scores)}).
               Generates final phenotypes, recorded in [this file]({path(input.mut_phenos_file)}).


            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"


rule get_UShER_tree:
    """Get UShER SARS-CoV-2 tree: https://github.com/yatisht/usher."""
    output:
        directory=directory(config['UShER_dir']),
        mat=config['UShER_tree'],
    shell:
        """
        wget \
            -r \
            -l 1 \
            -np \
            -nH \
            -R "index.html*" \
            -P {output.directory} \
            http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
        mv {output.directory}/goldenPath/wuhCor1/UShER_SARS-CoV-2/* {output.directory}
        rm -r {output.directory}/goldenPath
        """

rule get_UShER_refseq:
    """Get UShER reference sequence."""
    output: refseq=config['UShER_ref']
    shell:
        """
        efetch \
            -format gb \
            -db nuccore \
            -id NC_045512.2 \
            > {output.refseq}
        """

rule genomic_mutcounts:
    """Get counts of each genomic amino acid mutation."""
    input:
    	mat=config['UShER_tree']
    output:
    	genomic_mutcounts='results/genomic_mutcounts.csv'
    shell:
        # https://usher-wiki.readthedocs.io/en/latest/matUtils.html#summary
        """
        matUtils summary \
            -T 4 \
            -i {input.mat} \
            -m {output.genomic_mutcounts}
        """

rule collapse_scores:
    input:
        config['Titeseq_Kds_file'],
        config['expression_sortseq_file']
    output:
        config['final_variant_scores_mut_file'],
        md='results/summary/collapse_scores.md',
        md_files=directory('results/summary/collapse_scores_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='collapse_scores.Rmd',
        md='collapse_scores.md',
        md_files='collapse_scores_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_titrations:
    input:
        config['codon_variant_table_file_Wuhan_Hu_1'],
        config['codon_variant_table_file_E484K'],
        config['codon_variant_table_file_N501Y'],
        config['codon_variant_table_file_B1351'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file'],
        md='results/summary/compute_binding_Kd.md',
        md_files=directory('results/summary/compute_binding_Kd_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd.Rmd',
        md='compute_binding_Kd.md',
        md_files='compute_binding_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule calculate_expression:
    input:
        config['codon_variant_table_file_Wuhan_Hu_1'],
        config['codon_variant_table_file_E484K'],
        config['codon_variant_table_file_N501Y'],
        config['codon_variant_table_file_B1351'],
        config['variant_counts_file']
    output:
        config['expression_sortseq_file'],
        md='results/summary/compute_expression_meanF.md',
        md_files=directory('results/summary/compute_expression_meanF_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_expression_meanF.Rmd',
        md='compute_expression_meanF.md',
        md_files='compute_expression_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['codon_variant_table_file_Wuhan_Hu_1'],
        config['codon_variant_table_file_E484K'],
        config['codon_variant_table_file_N501Y'],
        config['codon_variant_table_file_B1351'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule get_mut_bind_expr:
    """Download SARS-CoV-2 mutation ACE2-binding and expression from URL."""
    output:
        file=config['mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind_expr_url'], output.file)
        

rule process_ccs_Wuhan_Hu_1:
    """Process the PacBio CCSs for Wuhan_Hu_1 background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_Wuhan_Hu_1'],
    	config['codon_variant_table_file_Wuhan_Hu_1'],
        nb_markdown=nb_markdown('process_ccs_Wuhan_Hu_1.ipynb')
    params:
        nb='process_ccs_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule process_ccs_E484K:
    """Process the PacBio CCSs for E484K background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_E484K'],
    	config['codon_variant_table_file_E484K'],
        nb_markdown=nb_markdown('process_ccs_E484K.ipynb')
    params:
        nb='process_ccs_E484K.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule process_ccs_N501Y:
    """Process the PacBio CCSs for N501Y background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_N501Y'],
    	config['codon_variant_table_file_N501Y'],
        nb_markdown=nb_markdown('process_ccs_N501Y.ipynb')
    params:
        nb='process_ccs_N501Y.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule process_ccs_B1351:
    """Process the PacBio CCSs for B1351 background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_B1351'],
    	config['codon_variant_table_file_B1351'],
        nb_markdown=nb_markdown('process_ccs_B1351.ipynb')
    params:
        nb='process_ccs_B1351.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

if config['seqdata_source'] == 'HutchServer':

    rule get_ccs:
        """Symbolically link CCS files."""
        input:
            ccs_fastq=lambda wildcards: (pacbio_runs
                                        .set_index('pacbioRun')
                                        .at[wildcards.pacbioRun, 'ccs']
                                        )
        output:
            ccs_fastq=os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz")
        run:
            os.symlink(input.ccs_fastq, output.ccs_fastq)

elif config['seqdata_source'] == 'SRA':
    raise RuntimeError('getting sequence data from SRA not yet implemented')

else:
    raise ValueError(f"invalid `seqdata_source` {config['seqdata_source']}")
