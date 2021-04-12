container: "docker://mfansler/scutr-quant:0.1.3"
configfile: "config.yaml"

import pandas as pd
import os
from glob import glob

# make sure the tmp directory exists
os.makedirs(config['tmp_dir'], exist_ok=True)

print("Loading sample data...")
samples = pd.read_csv(config['sample_file'], index_col='sample_id')
print("Loaded %d samples." % len(samples.index))

wildcard_constraints:
    sample_id=config['sample_regex']

rule all:
    input:
        config['final_output_file'],
        expand("qc/umi_count/{sample_id}.umi_count.html",
               sample_id=samples.index.values)


def get_file_type(wildcards):
    return "--bam" if samples.file_type[wildcards.sample_id] == 'bam' else ""

def get_sequence_files(wildcards):
    return samples.files[wildcards.sample_id].split(';')

rule kallisto_bus:
    input:
        idx=config['utrome_kdx'],
        files=get_sequence_files
    output:
        bus=temp("data/kallisto/{sample_id}/output.bus"),
        ec="data/kallisto/{sample_id}/matrix.ec",
        tx="data/kallisto/{sample_id}/transcripts.txt"
    params:
        tech=config['tech'],
        strand=config['strand'],
        bam=get_file_type
    threads: 16
    resources:
        mem_mb=1000
    shell:
        """
        outDir=$(dirname {output.bus})
        rm -r $outDir
        kallisto bus {params.bam} -t {threads} -i {input.idx} \\
        -x {params.tech} {params.strand} -o $outDir --verbose {input.files}
        """

rule bustools_sort:
    input:
        "data/kallisto/{sample_id}/output.bus"
    output:
        "data/kallisto/{sample_id}/output.sorted.bus"
    params:
        tmpDir=lambda wcs: config['tmp_dir'] + "/bs-utrome-sort" + wcs.sample_id
    threads: 4
    resources:
        mem_mb=2000
    shell:
        """
        bustools sort -t{threads} -T {params.tmpDir} -o {output} {input}
        """

rule bustools_whitelist:
    input:
        "data/kallisto/{sample_id}/output.sorted.bus"
    output:
        "data/kallisto/{sample_id}/whitelist.txt"
    shell:
        """
        bustools whitelist -o {output} {input}
        """

def get_whitelist(wildcards):
    if config['bx_whitelist'] == "":
        return "data/kallisto/%s/whitelist.txt" % wildcards.sample_id
    else:
        return config['bx_whitelist']
    
rule bustools_correct:
    input:
        bus="data/kallisto/{sample_id}/output.sorted.bus",
        bxs=get_whitelist
    output:
        temp("data/kallisto/{sample_id}/output.corrected.bus")
    shell:
        """
        bustools correct -w {input.bxs} -o {output} {input.bus}
        """

rule bustools_correct_sort:
    input:
        "data/kallisto/{sample_id}/output.corrected.bus"
    output:
        temp("data/kallisto/{sample_id}/output.corrected.sorted.bus")
    params:
        tmpDir=lambda wcs: config['tmp_dir'] + "/bs-utrome-sort" + wcs.sample_id
    threads: 4
    resources:
        mem_mb=2000
    shell:
        """
        bustools sort -t{threads} -T {params.tmpDir} -o {output} {input}
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

def get_input_busfile(wildcards):
    if config['correct_bus']:
        return "data/kallisto/%s/output.corrected.sorted.bus" % wildcards.sample_id
    else:
        return "data/kallisto/%s/output.sorted.bus" % wildcards.sample_id

rule bustools_count_txs:
    input:
        bus=get_input_busfile,
        txs="data/kallisto/{sample_id}/transcripts.txt",
        ec="data/kallisto/{sample_id}/matrix.ec",
        merge="data/utrs/utrome_tx_merge.tsv"
    output:
        mtx="data/kallisto/{sample_id}/utrome.txs.mtx",
        txs="data/kallisto/{sample_id}/utrome.txs.genes.txt",
        bxs="data/kallisto/{sample_id}/utrome.txs.barcodes.txt"
    shell:
        """
        filename={output.mtx}
        prefix=${{filename%.*}}
        bustools count -e {input.ec} -t {input.txs} -g {input.merge} -o $prefix --em --genecounts {input.bus}
        """

rule bustools_count_genes:
    input:
        bus=get_input_busfile,
        txs="data/kallisto/{sample_id}/transcripts.txt",
        ec="data/kallisto/{sample_id}/matrix.ec",
        merge="data/utrs/utrome_gene_merge.tsv"
    output:
        mtx="data/kallisto/{sample_id}/utrome.genes.mtx",
        txs="data/kallisto/{sample_id}/utrome.genes.genes.txt",
        bxs="data/kallisto/{sample_id}/utrome.genes.barcodes.txt"
    shell:
        """
        filename={output.mtx}
        prefix=${{filename%.*}}
        bustools count -e {input.ec} -t {input.txs} -g {input.merge} -o $prefix --em --genecounts {input.bus}
        """

rule report_umis_per_cell:
    input:
        bxs="data/kallisto/{sample_id}/utrome.txs.barcodes.txt",
        txs="data/kallisto/{sample_id}/utrome.txs.genes.txt",
        mtx="data/kallisto/{sample_id}/utrome.txs.mtx"
    params:
        min_umis=config['min_umis']
    output:
        "qc/umi_count/{sample_id}.umi_count.html"
    script:
        "scripts/report_umi_counts_per_cell.Rmd"

rule mtxs_to_sce:
    input:
        bxs=expand("data/kallisto/{sample_id}/utrome.txs.barcodes.txt", sample_id=samples.index.values),
        txs=expand("data/kallisto/{sample_id}/utrome.txs.genes.txt", sample_id=samples.index.values),
        mtxs=expand("data/kallisto/{sample_id}/utrome.txs.mtx", sample_id=samples.index.values),
        gtf=config['utrome_gtf']
    output:
        sce=config['final_output_file']
    params:
        genome=config['genome'],
        sample_ids=samples.index.values,
        min_umis=config['min_umis'],
        annots=config['annotation_file']
    resources:
        mem_mb=16000
    script:
        "scripts/mtxs_to_sce.R"

rule sce_txs_to_genes:
    input:
        sce="data/sce/{dataset}.txs.Rds"
    output:
        sce="data/sce/{dataset}.genes.Rds"
    resources:
        mem_mb=8000
    script:
        "scripts/sce_txs_to_genes.R"
