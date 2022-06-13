container: "docker://mfansler/scutr-quant:0.1.5"
configfile: "config.yaml"

from snakemake.io import load_configfile
import pandas as pd
import os
from sys import stderr

# print to stderr
def message(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

# make sure the tmp directory exists
os.makedirs(config['tmp_dir'], exist_ok=True)

# load target configuration
targets = load_configfile(config['targets_config'])

# check for valid target(s)
target_list = [config['target']] if isinstance(config['target'], str) else config['target']
if not set(target_list).issubset(targets.keys()):
    message("[Error] Some target(s) '%s' not found!" % config['target'])
    message("Available targets:", *targets.keys(), sep="\n  ")
    exit(1)

def get_target_file(key):
    def f(wcs):
        if not targets[wcs.target][key]:
            return []
        else:
            return targets[wcs.target]['path'] + targets[wcs.target][key]
    return f

# load sample data
message("[INFO] Loading sample data...")
samples = pd.read_csv(config['sample_file'], index_col='sample_id')
message("[INFO] Loaded %d samples." % len(samples.index))

wildcard_constraints:
    sample_id=config['sample_regex']

    
rule all:
    input:
        expand("data/sce/{target}/{dataset}.{output_type}.Rds",
               target=target_list,
               dataset=config['dataset_name'],
               output_type=config['output_type']),
        expand("qc/umi_count/{target}/{sample_id}.umi_count.html",
               target=target_list,
               sample_id=samples.index.values)

def get_file_type(wildcards):
    return "--bam" if samples.file_type[wildcards.sample_id] == 'bam' else ""

def get_sequence_files(wildcards):
    return samples.files[wildcards.sample_id].split(';')

rule kallisto_bus:
    input:
        idx=get_target_file('kdx'),
        files=get_sequence_files
    output:
        bus=temp("data/kallisto/{target}/{sample_id}/output.bus"),
        ec="data/kallisto/{target}/{sample_id}/matrix.ec",
        tx="data/kallisto/{target}/{sample_id}/transcripts.txt"
    params:
        tech=config['tech'],
        strand=config['strand'],
        bam=get_file_type
    threads: 16
    resources:
        mem_mb=1000
    conda: "envs/kallisto-bustools.yaml"
    shell:
        """
        outDir=$(dirname {output.bus})
        rm -r $outDir
        kallisto bus {params.bam} -t {threads} -i {input.idx} \\
        -x {params.tech} {params.strand} -o $outDir --verbose {input.files}
        """

rule bustools_sort:
    input:
        "data/kallisto/{target}/{sample_id}/output.bus"
    output:
        "data/kallisto/{target}/{sample_id}/output.sorted.bus"
    params:
        tmpDir=lambda wcs: config['tmp_dir'] + "/bs-utrome-sort" + wcs.sample_id
    threads: 4
    resources:
        mem_mb=2000
    conda: "envs/kallisto-bustools.yaml"
    shell:
        """
        bustools sort -t{threads} -T {params.tmpDir} -o {output} {input}
        """

rule bustools_whitelist:
    input:
        "data/kallisto/{target}/{sample_id}/output.sorted.bus"
    output:
        "data/kallisto/{target}/{sample_id}/whitelist.txt"
    conda: "envs/kallisto-bustools.yaml"
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
        bus="data/kallisto/{target}/{sample_id}/output.sorted.bus",
        bxs=get_whitelist
    output:
        temp("data/kallisto/{target}/{sample_id}/output.corrected.bus")
    conda: "envs/kallisto-bustools.yaml"
    shell:
        """
        bustools correct -w {input.bxs} -o {output} {input.bus}
        """

rule bustools_correct_sort:
    input:
        "data/kallisto/{target}/{sample_id}/output.corrected.bus"
    output:
        temp("data/kallisto/{target}/{sample_id}/output.corrected.sorted.bus")
    params:
        tmpDir=lambda wcs: config['tmp_dir'] + "/bs-utrome-sort" + wcs.target + "-" + wcs.sample_id
    threads: 4
    resources:
        mem_mb=2000
    conda: "envs/kallisto-bustools.yaml"
    shell:
        """
        bustools sort -t{threads} -T {params.tmpDir} -o {output} {input}
        """

rule generate_tx_merge:
    input:
        tsv=get_target_file('merge_tsv')
    output:
        "data/utrs/{target}/tx_merge.tsv"
    shell:
        """
        tail -n+2 {input.tsv} | cut -f1,2 > {output}
        """

rule generate_gene_merge:
    input:
        tsv=get_target_file('merge_tsv')
    output:
        "data/utrs/{target}/gene_merge.tsv"
    shell:
        """
        tail -n+2 {input.tsv} | cut -f1,3 > {output}
        """

def get_input_busfile(wildcards):
    if config['correct_bus']:
        return "data/kallisto/%s/%s/output.corrected.sorted.bus" % (wildcards.target, wildcards.sample_id)
    else:
        return "data/kallisto/%s/%s/output.sorted.bus" % (wildcards.target, wildcards.sample_id)

rule bustools_count_txs:
    input:
        bus=get_input_busfile,
        txs="data/kallisto/{target}/{sample_id}/transcripts.txt",
        ec="data/kallisto/{target}/{sample_id}/matrix.ec",
        merge="data/utrs/{target}/tx_merge.tsv"
    output:
        mtx="data/kallisto/{target}/{sample_id}/txs.mtx",
        txs="data/kallisto/{target}/{sample_id}/txs.genes.txt",
        bxs="data/kallisto/{target}/{sample_id}/txs.barcodes.txt"
    conda: "envs/kallisto-bustools.yaml"
    shell:
        """
        filename={output.mtx}
        prefix=${{filename%.*}}
        bustools count -e {input.ec} -t {input.txs} -g {input.merge} -o $prefix --em --genecounts {input.bus}
        """

rule bustools_count_genes:
    input:
        bus=get_input_busfile,
        txs="data/kallisto/{target}/{sample_id}/transcripts.txt",
        ec="data/kallisto/{target}/{sample_id}/matrix.ec",
        merge="data/utrs/{target}/gene_merge.tsv"
    output:
        mtx="data/kallisto/{target}/{sample_id}/genes.mtx",
        txs="data/kallisto/{target}/{sample_id}/genes.genes.txt",
        bxs="data/kallisto/{target}/{sample_id}/genes.barcodes.txt"
    conda: "envs/kallisto-bustools.yaml"
    shell:
        """
        filename={output.mtx}
        prefix=${{filename%.*}}
        bustools count -e {input.ec} -t {input.txs} -g {input.merge} -o $prefix --em --genecounts {input.bus}
        """

rule report_umis_per_cell:
    input:
        bxs="data/kallisto/{target}/{sample_id}/txs.barcodes.txt",
        txs="data/kallisto/{target}/{sample_id}/txs.genes.txt",
        mtx="data/kallisto/{target}/{sample_id}/txs.mtx"
    params:
        min_umis=config['min_umis']
    output:
        "qc/umi_count/{target}/{sample_id}.umi_count.html"
    conda: "envs/rmd-reporting.yaml"
    script:
        "scripts/report_umi_counts_per_cell.Rmd"

rule mtxs_to_sce_txs:
    input:
        bxs=expand("data/kallisto/{target}/{sample_id}/txs.barcodes.txt",
                   sample_id=samples.index.values, allow_missing=True),
        txs=expand("data/kallisto/{target}/{sample_id}/txs.genes.txt",
                   sample_id=samples.index.values, allow_missing=True),
        mtxs=expand("data/kallisto/{target}/{sample_id}/txs.mtx",
                    sample_id=samples.index.values, allow_missing=True),
        gtf=get_target_file('gtf'),
        tx_annots=get_target_file('tx_annots')
    output:
        sce="data/sce/{target}/%s.txs.Rds" % config['dataset_name']
    params:
        genome=lambda wcs: targets[wcs.target]['genome'],
        sample_ids=samples.index.values,
        min_umis=config['min_umis'],
        cell_annots=config['cell_annots']
    resources:
        mem_mb=16000
    conda: "envs/bioconductor-sce.yaml"
    script:
        "scripts/mtxs_to_sce_txs.R"

rule mtxs_to_sce_genes:
    input:
        bxs=expand("data/kallisto/{target}/{sample_id}/genes.barcodes.txt",
                   sample_id=samples.index.values, allow_missing=True),
        genes=expand("data/kallisto/{target}/{sample_id}/genes.genes.txt",
                     sample_id=samples.index.values, allow_missing=True),
        mtxs=expand("data/kallisto/{target}/{sample_id}/genes.mtx",
                    sample_id=samples.index.values, allow_missing=True),
        gtf=get_target_file('gtf'),
        gene_annots=get_target_file('gene_annots')
    output:
        sce="data/sce/{target}/%s.genes.Rds" % config['dataset_name']
    params:
        genome=lambda wcs: targets[wcs.target]['genome'],
        sample_ids=samples.index.values,
        min_umis=config['min_umis'],
        cell_annots=config['cell_annots']
    resources:
        mem_mb=16000
    conda: "envs/bioconductor-sce.yaml"
    script:
        "scripts/mtxs_to_sce_genes.R"
