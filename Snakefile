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

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        get_mut_bind_expr=config['mut_bind_expr'],
        get_delta_mut_bind_expr=config['delta_mut_bind_expr'],
        get_mut_antibody_escape=config['mut_antibody_escape'],
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
        sars2_subs=config['UShER_annotated_subs'],
        parsed_subs_N501Y=config['UShER_parsed_subs_N501Y'],
        epistatic_shifts='results/summary/epistatic_shifts.md'
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

            1. Get prior Wuhan-1 RBD DMS mutation-level [binding and expression data]({path(input.get_mut_bind_expr)}) and Delta VOC RBD DMS mutation-level [binding and expression data]{path(input.get_delta_mut_bind_expr)}.
            
            2. Process PacBio CCSs for each background: [Wuhan_Hu_1]({path(input.process_ccs_Wuhan_Hu_1)}), [E484K]({path(input.process_ccs_E484K)}), [N501Y]({path(input.process_ccs_N501Y)}), [B.1.351]({path(input.process_ccs_B1351)}). Creates barcode-variant lookup tables for each background: [Wuhan_Hu_1]({path(input.barcode_variant_table_Wuhan_Hu_1)}), [E484K]({path(input.barcode_variant_table_E484K)}), [N501Y]({path(input.barcode_variant_table_N501Y)}), [B.1.351]({path(input.barcode_variant_table_B1351)}).
            
            3. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            4. [Fit titration curves]({path(input.fit_titrations)}) to calculate per-barcode K<sub>D</sub>, recorded in [this file]({path(input.variant_Kds_file)}).
            
            5. [Analyze Sort-seq]({path(input.calculate_expression)}) to calculate per-barcode RBD expression, recorded in [this file]({path(input.variant_expression_file)}).
            
            6. [Derive final genotype-level phenotypes from replicate barcoded sequences]({path(input.collapse_scores)}).
               Generates final phenotypes, recorded in [this file]({path(input.mut_phenos_file)}).
               
            7. Download trees and reference files for analysis of the mutation-annotated tree provided in [UShER](https://github.com/yatisht/usher), and parse substitution occurrence counts by background: [N501Y]({path(input.parsed_subs_N501Y)})
            
            8. [Analyze patterns of epistasis in the DMS data and in SARS-CoV-2 genomic data]({path(input.epistatic_shifts)}).


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

rule epistatic_shifts:
    input:
        config['final_variant_scores_mut_file'],
        config['mut_antibody_escape'],
        config['UShER_parsed_subs_N501Y']
    output:
        config['JSD_v_WH1_file'],
        md='results/summary/epistatic_shifts.md',
        md_files=directory('results/summary/epistatic_shifts_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='epistatic_shifts.Rmd',
        md='epistatic_shifts.md',
        md_files='epistatic_shifts_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule get_UShER_tree:
    """Get UShER SARS-CoV-2 tree: https://github.com/yatisht/usher."""
    output:
        mat=config['UShER_tree'],
        nwk=config['UShER_tree_nwk']
    params:
        directory=directory(config['UShER_dir']),
        nwk_gz=config['UShER_tree_nwk_gz']
    shell:
        """
        mkdir -p {params.directory}
        wget \
            -r \
            -l 1 \
            -np \
            -nH \
            -R "index.html*" \
            -P {params.directory} \
            http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
        mv {params.directory}/goldenPath/wuhCor1/UShER_SARS-CoV-2/* {params.directory}
        rm -r {params.directory}/goldenPath
        gunzip {params.nwk_gz}
        """

rule get_UShER_gtf_file:
    """Get UShER reference gene annotations."""
    output:
    	gtf=config['UShER_gtf']
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/genes/ncbiGenes.gtf.gz
        gunzip ncbiGenes.gtf.gz
        mv ./ncbiGenes.gtf {output.gtf}
        """

rule get_UShER_refseq_file:
    """Get UShER reference sequence fasta."""
    output:
    	refseq=config['UShER_ref'],
    run:
        urllib.request.urlretrieve("https://raw.githubusercontent.com/yatisht/usher/master/test/NC_045512v2.fa", output.refseq)
    

rule enumerate_UShER_subs:
    """Get counts of each amino acid mutation on the UShER MAT."""
    input:
    	mat=config['UShER_tree'],
    	gtf=config['UShER_gtf'],
    	refseq=config['UShER_ref']
    output:
    	subs=config['UShER_annotated_subs']
    shell:
        # https://usher-wiki.readthedocs.io/en/latest/matUtils.html#summary
        """
        matUtils summary --translate {output.subs} -i {input.mat} -g {input.gtf} -f {input.refseq}
        """

rule parse_UShER_subs_N501Y:
    """Get counts of substitutions on N501 versus Y501 backgrounds."""
    input:
    	nwk=config['UShER_tree_nwk'],
    	subs=config['UShER_annotated_subs']
    output:
    	parsed_subs_N501Y=config['UShER_parsed_subs_N501Y']
    shell:
        """
        python ./scripts/Will_search_mat_211123.py -t {input.nwk} -a {input.subs} -s N501Y --skip_transitions=False -o {output.parsed_subs_N501Y}
        """


rule collapse_scores:
    input:
        config['Titeseq_Kds_file'],
        config['expression_sortseq_file'],
        config['delta_mut_bind_expr']
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

rule get_delta_mut_bind_expr:
    """Download SARS-CoV-2 Delta VOC mutation-level ACE2-binding and expression data."""
    output:
        file=config['delta_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['delta_mut_bind_expr_url'], output.file)

rule get_mut_antibody_escape:
    """Download SARS-CoV-2 mutation antibody-escape data from URL."""
    output:
        file=config['mut_antibody_escape']
    run:
        urllib.request.urlretrieve(config['mut_antibody_escape_url'], output.file)
        
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
