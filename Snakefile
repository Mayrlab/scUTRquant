configfile: "config.yaml"

import pandas as pd
import os
from glob import glob

# make sure the tmp directory exists
os.makedirs(config['tmp_dir'], exist_ok=True)

print("Loading sample data...")
samples = pd.read_csv(config['sample_file'], sep='\t', index_col='sample_id')
print("Loaded %d samples." % length(samples.index))

wildcard_constraints:
    srr="SRR\d+",
    sample_id=config['sample_regex']

rule all:
    input:
        config['final_output_file'],
        #expand("qc/umi_count/{sample_id}.UMI_Count.html",
        #       sample_id=samples.index.values)


def get_file_type(wildcards):
    return "--bam" if samples.file_type.values[wildcards.sample_id] == 'bam' else ""

def get_sequence_files(wildcards):
    return samples.files.values[wildcars.sample_id].split(';')

rule kallisto_bus:
    input:
        idx=config['utrome_kdx'],
        files=get_sequence_files
    output:
        bus="data/kallisto/{sample_id}/output.bus",
        ec="data/kallisto/{sample_id}/matrix.ec",
        tx="data/kallisto/{sample_id}/transcripts.txt"
    params:
        tech=config['tech'],
        strand=config['strand'],
        bam=get_file_type
    threads: 16
    resources:
        mem=1
    shell:
        """
        outDir=$(dirname {output.bus})
        rm -r $outDir
        kallisto bus {params.bam} -t {threads} -i {input.idx} \\
        -x {params.tech} {params.strand} -o $outDir --verbose {input.files}
        """

rule bustools_correct:
    input:
        bus="data/kallisto/{sample_id}/output.bus",
        bxs=config['bx_whitelist']
    output:
        temp("data/kallisto/{sample_id}/output.corrected.bus")
    shell:
        """
        bustools correct -w {input.bxs} -o {output} {input.bus}
        """

rule bustools_sort:
    input:
        "data/kallisto/{sample_id}/output.corrected.bus"
    output:
        "data/kallisto/{sample_id}/output.corrected.sorted.bus"
    params:
        tmpDir=lambda wcs: config['tmp_dir'] + "/bs-utrome-sort" + wcs.sample_id
    threads: 8
    resources:
        mem=2
    shell:
        """
        bustools sort -t{threads} -m14G -T {params.tmpDir} -o {output} {input}
        """

rule generate_tx_merge:
    input:
        config['utrome_merge']
    output:
        "data/utrs/utrome_tx_merge.tsv"
    shell:
        """
        tail -n+2 {input} | cut -f1,2 > {output}
        """

rule generate_gene_merge:
    input:
        config['utrome_merge']
    output:
        "data/utrs/utrome_gene_merge.tsv"
    shell:
        """
        tail -n+2 {input} | cut -f1,3 > {output}
        """

rule bustools_count_txs:
    input:
        bus="data/kallisto/{sample_id}/output.corrected.sorted.bus",
        txs="data/kallisto/{sample_id}/transcripts.txt",
        ec="data/kallisto/{sample_id}/matrix.ec",
        merge="data/utrs/utrome_tx_merge.tsv"
    output:
        "data/kallisto/{sample_id}/utrome.txs.mtx"
    shell:
        """
        filename={output}
        prefix=${{filename%.*}}
        bustools count -e {input.ec} -t {input.txs} -g {input.merge} -o $prefix --em --genecounts {input.bus}
        """

rule bustools_count_genes:
    input:
        bus="data/kallisto/{sample_id}/output.corrected.sorted.bus",
        txs="data/kallisto/{sample_id}/transcripts.txt",
        ec="data/kallisto/{sample_id}/matrix.ec",
        merge="data/utrs/utrome_gene_merge.tsv"
    output:
        "data/kallisto/{sample_id}/utrome.genes.mtx"
    conda:
        "envs/bustools.yaml"
    shell:
        """
        filename={output}
        prefix=${{filename%.*}}
        bustools count -e {input.ec} -t {input.txs} -g {input.merge} -o $prefix --em --genecounts {input.bus}
        """

rule report_umis_per_cell:
    input:
        bxs="data/kallisto/{sample_id}/utrome.txs.barcodes.txt",
        txs="data/kallisto/{sample_id}/utrome.txs.genes.txt",
        mtx="data/kallisto/{sample_id}/utrome.txs.mtx",
        annots=config['annotation_file']
    output:
        "qc/umi_count/{sample_id}.umi_count.html"
    conda: "envs/sce.yaml"
    script:
        "scripts/report_umi_counts_per_cell.Rmd"

rule mtxs_to_sce:
    input:
        bxs=expand("data/kallisto/{sample_id}/utrome.txs.barcodes.txt", sample_id=samples.index.values),
        txs=expand("data/kallisto/{sample_id}/utrome.txs.genes.txt", sample_id=samples.index.values),
        mtxs=expand("data/kallisto/{sample_id}/utrome.txs.mtx", sample_id=samples.index.values),
        gtf=config['utrome_gtf'],
        annots=config['annotation_file']
    output:
        sce=config['final_output_file']
    params:
        genome=config['genome'],
        sample_ids=samples.index.values
    resources:
        mem=16
    script:
        "scripts/mtxs_to_sce.R"
