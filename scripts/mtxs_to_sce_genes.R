#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(plyranges)
library(Matrix)
library(SingleCellExperiment)

################################################################################
# Fake Data (Interactive Testing)
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(
            bxs=c("data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.barcodes.txt"),
            genes=c("data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.genes.txt"),
            mtxs=c("data/kallisto/utrome_hg38_v1/pbmc_1k_v2_fastq/genes.mtx"),
            gtf="extdata/targets/utrome_hg38_v1/utrome.e30.t5.gc39.pas3.f0.9999.w500.gtf",
            gene_annots=NULL,
            cell_annots="examples/pbmc_1k_v2_fastq/annots.csv"),
        output=list(sce="data/sce/utrome_hg38_1/test.utrome.genes.Rds"),
        params=list(genome='hg38',
                    sample_ids=c("pbmc_1k_v2_fastq"),
                    min_umis="500",
                    cell_annots_key="cell_id"))
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

################################################################################
# Load Data
################################################################################

gr_gtf <- read_gff(snakemake@input$gtf, genome_info=snakemake@params$genome,
                   col_names=c('type', 'transcript_id', 'gene_id',
                               'gene_name', 'transcript_name', 'utr_name', 'utr_rank'))
grl_genes <- gr_gtf %>%
    filter(type == 'transcript') %>%

    ## remove unneeded or missing columns
    select(-c('type', 'gene_name')) %>%
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
df_cell_annots <- snakemake@input$cell_annots %>% {
    if (is.null(.)) {
        tibble(!!(snakemake@params$cell_annots_key):=character(0))
    } else {
        read_csv(.)
    }
}

## Load MTXs Function
load_mtx_to_sce <- function (mtxFile, bxFile, geneFile, sample_id) {
    ## TODO: Add annotation-based filtering
    ## valid_cells <- filter(df_annots, sample_id == sample_id) %$% cell_id
    readMM(mtxFile) %>%
        as("CsparseMatrix") %>%
        `dimnames<-`(list(
            cell_id=str_c(sample_id, "_", read_lines(bxFile)),
            gene_id=read_lines(geneFile))) %>%
        t %>%
        { .[, colSums(.) >= as.integer(snakemake@params$min_umis)] } %>%
        { SingleCellExperiment(assays=list(counts=.),
                               colData=DataFrame(cell_id=colnames(.),
                                                 sample_id=sample_id,
                                                 row.names=colnames(.))) }
}


sce <- mapply(load_mtx_to_sce,
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

saveRDS(sce, snakemake@output$sce)
