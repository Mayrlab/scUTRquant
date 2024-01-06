#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(plyranges)
library(Matrix)
library(HDF5Array)
library(SingleCellExperiment)
library(BiocParallel)

################################################################################
# Fake Data (Interactive Testing)
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               params='list', threads='numeric'))
    snakemake <- Snakemake(
        input=list(
            bxs=c("data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.barcodes.txt",
                  "data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.barcodes.txt"),
            genes=c("data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.genes.txt",
                  "data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.genes.txt"),
            mtxs=c("data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.mtx",
                   "data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.mtx"),
            gtf="extdata/targets/utrome_hg38_v1/utrome.e30.t5.gc39.pas3.f0.9999.w500.gtf",
            gene_annots=NULL,
            cell_annots="examples/pbmc_1k_v2_fastq/annots.csv"),
        output=list(sce="data/sce/utrome_hg38_v1/test.utrome.genes.se.rds",
                    hdf5="data/sce/utrome_hg38_v1/test.utrome.genes.assays.h5"),
        params=list(genome='hg38',
                    sample_ids=c("pbmc_1k_v2_fastq_1", "pbmc_1k_v2_fastq_2"),
                    min_umis="500",
                    cell_annots_key="cell_id",
                    exclude_unannotated_cells=FALSE,
                    tmp_dir=tempdir(check=TRUE),
                    use_hdf5=TRUE),
        threads=4L)
}

################################################################################
# Helper Functions
################################################################################
is_na_all <- function (v) { all(is.na(v)) }

drop_na_mcols <- function (gr) {
    cols_to_rm <- mcols(gr) %>%
        apply(2, is_na_all) %>%
        which %>% names
    select(gr, -cols_to_rm)
}

filter_min_umis <- function (x) {
    x[, colSums(x) >= as.integer(snakemake@params$min_umis)]
}

optional_hdf5_convert <- if (snakemake@params$use_hdf5) {
    dir.create(snakemake@params$tmp_dir, recursive=TRUE)
    setHDF5DumpDir(snakemake@params$tmp_dir)
    writeTENxMatrix
} else { identity }

################################################################################
# Load Data
################################################################################

gr_gtf <- read_gff(snakemake@input$gtf, genome_info=snakemake@params$genome,
                   col_names=c('type', 'transcript_id', 'gene_id',
                               'gene_name', 'transcript_name', 'utr_name', 'utr_rank'))

grl_genes <- gr_gtf %>%
    filter(type == 'transcript') %>%

    ## remove unneeded or missing columns
    select(-all_of(c('type', 'gene_name'))) %>%
    drop_na_mcols %>%

    ## add rownames
    `names<-`(.$transcript_id) %>%

    ## group by gene_id
    split(.$gene_id)

## Row Annotations (UTRs)
has_row_annots <- !is.null(snakemake@input$gene_annots)
if (has_row_annots) {
    df_gene_annots <- readRDS(snakemake@input$gene_annots)
} else {
    df_gene_annots <- gr_gtf %>%

        ## keep only transcripts
        filter(type == 'gene') %>%

        ## remove unneeded or missing columns
        select(-c('type', 'transcript_id', 'transcript_name', 'utr_name', 'utr_rank')) %>%
        drop_na_mcols %>%

        ## add rownames
        `names<-`(.$gene_id) %>%

        mcols() %>%

        unique()
}

## Column Annotations
has_cell_annots <- !is.null(snakemake@input$cell_annots)
df_cell_annots <- if (has_cell_annots) {
    read_csv(snakemake@input$cell_annots)
} else {
    tibble(!!(snakemake@params$cell_annots_key):=character(0))
}

optional_filter_unannotated <- if (has_cell_annots && snakemake@params$exclude_unannotated_cells) {
    function (x) {
        x[, colnames(x) %in% df_cell_annots[[snakemake@params$cell_annots_key]]]
    }
} else { identity }

## Load MTXs Function
load_mtx_to_sce <- function (mtxFile, bxFile, geneFile, sample_id) {
    readMM(mtxFile) %>%
        as("CsparseMatrix") %>%
        `dimnames<-`(list(
            cell_id=str_c(sample_id, "_", read_lines(bxFile)),
            gene_id=read_lines(geneFile))) %>%
        t %>%
        optional_filter_unannotated %>%
        filter_min_umis %>%
        optional_hdf5_convert %>%
        { SingleCellExperiment(assays=list(counts=.),
                               colData=DataFrame(cell_id=colnames(.),
                                                 sample_id.sq=sample_id,
                                                 row.names=colnames(.))) }
}

register(MulticoreParam(snakemake@threads))

sce <- bpmapply(load_mtx_to_sce,
                mtxFile=snakemake@input$mtxs,
                bxFile=snakemake@input$bxs,
                geneFile=snakemake@input$genes,
                sample_id=snakemake@params$sample_ids) %>%
    do.call(what=cbind)

################################################################################
# Attach Annotations
################################################################################

colData(sce) %<>%
    as_tibble %>%
    left_join(df_cell_annots, by=c("cell_id"=snakemake@params$cell_annots_key)) %>%
    set_rownames(.$cell_id) %>%
    DataFrame %>%
    `[`(colnames(sce),)

## attach location info
rowRanges(sce) <- grl_genes[rownames(sce), ]

rowData(sce) <- df_gene_annots[rownames(sce), ]

################################################################################
# Export SCE
################################################################################

if (snakemake@params$use_hdf5) {
    chunkdim <- c(dim(sce)[[1]], min(1e3L, dim(sce)[[2]]))
    saveHDF5SummarizedExperiment(sce, dir=dirname(snakemake@output$sce),
                                 prefix=str_remove(basename(snakemake@output$sce), "se.rds$"),
                                 chunkdim=chunkdim, as.sparse=TRUE, verbose=TRUE)
    unlink(snakemake@params$tmp_dir, recursive=TRUE, force=TRUE)
} else {
    saveRDS(sce, snakemake@output$sce)
}
