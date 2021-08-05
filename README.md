# scUTRquant
A bioinformatics pipeline for single-cell 3' UTR isoform quantification.

# Overview
The **scUTRquant** pipeline builds on `kallisto bus` to provide a reusable tool for 
quantifying 3' UTR isoforms from 3'-end tag-based scRNA-seq datasets. The pipeline
is based on Snakemake and reproducibily bundles all processing software in a single
Docker image to abstract away software installations. It includes a prebuilt 
reference UTRome for **mm10**, which ensures a consistent set of features (3' UTR isoforms)
across different runs. In total, this provides a rapid pipeline for recovering
3' UTR isoform counts from common scRNA-seq datasets.

## Inputs
The pipeline takes as input:

- set of FASTQ or BAM (CellRanger output) files from scRNA-seq experiments
- kallisto index of UTRome (**mm10** provided)
- GTF annotation of UTRome (**mm10** provided)
- TSV merge annotation (**mm10** provided)
- YAML configuration file that controls pipeline parameters
- CSV sample sheet detailing the FASTQ/BAM files to be processed
- barcode whitelist (optional)
- CSV of cell annotations (optional)

**Note on UTRome Index**

The pipeline includes code to download a prebuilt **mm10** UTRome GTF and kallisto index, 
as well as the 10X barcode whitelists (versions 1-3). This prebuilt index was generated
by augmenting the protein coding transcripts that have verified 3' ends in the GENCODE 
vM21 annotation with high-confidence cleavage sites called from the Mouse Cell Atlas 
dataset. This augmented transcriptome was then truncated to include only the last 
500 nts of each transcript and then deduplicated. Finally, the merge file contains
information on transcripts whose cleavage sites differ by fewer than 200 nts, which 
corresponds to the empirical resolution limit for `kallisto` quantification as 
determined by simulations.

## Outputs
The primary output of the pipeline is a Bioconductor `SingleCellExperiment` object.
The `counts` in the object is a sparse `Matrix` of 3' UTR isoform counts; the `rowRanges` 
is a `GenomicRanges` of the 3' UTR isoforms; the `rowData` is a `DataFrame` with additional
information about 3' UTR isoforms; and the `colData` is a `DataFrame` populated with sample 
metadata and optional user-provide cell annotations.

To assist users in quality control, the pipeline additionally generates HTML reports 
for each sample.

The pipeline is configured to retain intermediate files, such as BUS and MTX files.
Advanced users can readily customize the pipeline to only generate the files they 
require. For example, users who prefer to work with alternative scRNA-seq data structures,
such as those used in Scanpy or Seurat, may wish to terminate the pipeline at MTX 
generation.

# Setup
## Requirements
The following must be installed to execute the pipeline:

 - [Singularity](https://singularity.lbl.gov/index.html)
 - [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

This pipeline has only been tested on Linux. While it may be possible to run on
other systems, note that Snakemake is not available for Windows, and Singularity
only has native support for Linux.

## Installation
1. Clone the repository.
    ```
    git clone git@github.com:mfansler/scutr-quant.git
    ```

2. Download the UTRome annotation, kallisto index, and merge file.
    ```
    cd scutr-quant/extdata
    . download_utrome.sh
    ```
    **Reuse Tip:** For use across multiple projects, it is recommended to centralize 
    these files and change the entries in the `configfile` for `utrome_gtf`,
    `utrome_kdx`, and `utrome_merge` to point to the central location. In that
    case, one does not need to redownload the files.

3. (Optional) Download the barcode whitelists.
    ```
    cd scutr-quant/extdata/bxs
    . download_10X_whitelists.sh
    ```
    **Reuse Tip:** Similar to the UTRome files, these can also be centralized
    and referenced by the `bx_whitelist` variable in the `configfile`.

# Running Examples
Examples are provided in the `scutr-quant/examples` folder. Each includes a script
for downloading the raw data, a `sample_sheet.csv` formatted for use in the pipeline,
and a `config.yaml` file for running the pipeline.

Note that the `config.yaml` uses paths relative to the `scutr-quant` folder.

## 1K Neurons (10xv3) - BAM

1. Download the raw data.
    ```
    cd scutr-quant/examples/neuron_1k_v3_bam/
    . download.sh
    ```

2. Run the pipeline.
    ```
    cd scutr-quant
    snakemake --use-singularity --configfile examples/neuron_1k_v3_bam/config.yaml
    ```

## 1K Neurons (10xv3) - FASTQ

1. Download the raw data.
    ```
    cd scutr-quant/examples/neuron_1k_v3_fastq/
    . download.sh
    ```

2. Run the pipeline.
    ```
    cd scutr-quant
    snakemake --use-singularity --configfile examples/neuron_1k_v3_fastq/config.yaml
    ```

# File Specifications
## Configuration File

The Snakemake `configfile` specifies the parameters used to run the
pipeline. The following keys are expected:

 - `dataset_name`: name used in the final `SingleCellExperiment` object
 - `tmp_dir`: path to use for temporary files
 - `sample_file`: CSV-formatted file listing the samples to be processed
 - `utrome_gtf`: GTF annotation of UTRome; used in annotating rows
 - `utrome_kdx`: Kallisto index for UTRome
 - `utrome_merge`: TSV for merging features (isoforms)
 - `genome`: corresponds to genome that `utrome_gtf` references (only `mm10` currently supported)
 - `sample_regex`: regular expression used to match sample IDs; including a specific
     regex helps to constrain Snakemake's DAG-generation
 - `output_type`: a list of outputs, including `"txs"` and/or `"genes"`
 - `tech`: argument to `kallisto bus` indicating the scRNA-seq technology; see
     [the documentation](https://pachterlab.github.io/kallisto/manual#bus) for supported values
 - `strand`: argument to `kallisto bus` indicating the orientation of sequence reads
     with respect to transcripts; all 10X 3'-end libraries use `--fr-stranded`;
     omitting this argument eliminates the ability to correctly assign reads to
     transcripts when opposing stranded genes overlap
 - `bx_whitelist`: file of valid barcodes used in `bustools correct`
 - `min_umis`: minimum number of UMIs per cell; cells below this threshold are excluded
 
### Default Values

Snakemake can draw values for `config` in three ways:

 1. `scutr-quant/config.yaml`: This file is listed as the `configfile` in the Snakefile. 
 2. `--configfile config.yaml`: The file provided at the commandline.
 3. `--config argument=value`: A specific value for an argument 
 
This list runs from lowest to highest precedence. Configuration values that do not differ from those in `scutr-quant/config.yaml` can be left unspecfied in the YAML given by the `--configfile` argument. That is, one can use the `scutr-quant/config.yaml` to define shared settings, and only list dataset-specific config values in the dataset's YAML.

## Sample File

The `sample_file` provided in the Snakemake configuration is expected to be a CSV
with at least the following columns:

 - `sample_id`: a unique identifier for the sample; used in file names and paths
     of intermediate files derived from the sample
 - `file_type`: indicates whether sample input is `'bam'` or `'fastq'` format
 - `files`: a semicolon-separated list of files; for multi-run (e.g., multi-lane)
     samples, the files must have the order:
     ```
     lane1_R1.fastq;lane1_R2.fastq;lane2_R1.fastq;lane2_R2.fastq;...
     ```
 
# Customization

The rules in the `Snakefile` include `threads` and `resources` arguments per rule. These values are compatible for use with [Snakemake profiles](https://github.com/Snakemake-Profiles) for cluster deployment. The current defaults will attempt to use up to 16 threads and 16GB of memory in the `kallisto bus` step. Please adjust to fit the resources available on the deployment cluster. We strongly recommend that cluster profiles include both `--use-singularity` and `--use-conda` flags by default. Following this recommendation, an example run, for instance on **neuron_1k_v3_fastq**, with profile name `profile_name`, would take the form:

```bash
snakemake --profile profile_name --configfile examples/neuron_1k_v3_fastq/config.yaml
```
