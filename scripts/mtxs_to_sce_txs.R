#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(readr)
library(stringr)
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
            bxs=c("data/kallisto/day5_DP_1/utrome.txs.barcodes.txt",
                  "data/kallisto/day6_SP_1/utrome.txs.barcodes.txt"),
            txs=c("data/kallisto/day5_DP_1/utrome.txs.genes.txt",
                  "data/kallisto/day6_SP_1/utrome.txs.genes.txt"),
            mtxs=c("data/kallisto/day5_DP_1/utrome.txs.mtx",
                   "data/kallisto/day6_SP_1/utrome.txs.mtx"),
            gtf="extdata/gff/adult.utrome.e3.t200.f0.999.w500.gtf",
            tx_annots="extdata/atlas/utrome_txs_annotation.Rds"),
        output=list(sce="data/sce/test.utrome.txs.Rds"),
        params=list(genome='mm10',
                    sample_ids=c("day5_DP_1", "day6_SP_1"),
                    min_umis="500",
                    cell_annots="metadata/dahlin18_wolf19_annots.csv"))
}

################################################################################
                                        # Load Data
################################################################################

## Row Annotations (UTRs)
has_row_annots <- !is.null(snakemake@input$tx_annots)
if (has_row_annots) {
    df_tx_annots <- readRDS(snakemake@input$tx_annots)
}

## Row Ranges (
gr_txs <- read_gff(snakemake@input$gtf, genome_info=snakemake@params$genome,
                    col_names=c('type', 'transcript_id')) %>%

    ## keep only transcripts
    filter(type == 'transcript') %>%
    select(-c('type')) %>%

    ## add rownames
    `names<-`(.$transcript_id)

## Column Annotations
df_cell_annots <- snakemake@params$cell_annots %>% {
    if (is.null(.)) {
        tibble(cell_id=character(0))
    } else {
        read_csv(.) 
    }
}

## Load MTXs Function
load_mtx_to_sce <- function (mtxFile, bxFile, txFile, sample_id) {
    ## TODO: Add annotation-based filtering
    ## valid_cells <- filter(df_annots, sample_id == sample_id) %$% cell_id
    readMM(mtxFile) %>%
        as("CsparseMatrix") %>%
        `dimnames<-`(list(
            cell_id=str_c(sample_id, "_", read_lines(bxFile)),
            transcript_id=read_lines(txFile))) %>%
        t %>%
        { .[, colSums(.) >= as.integer(snakemake@params$min_umis)] } %>%
        { SingleCellExperiment(assays=list(counts=.),
                               colData=DataFrame(cell_id=colnames(.), sample_id=sample_id, row.names=colnames(.))) }
}


sce <- mapply(load_mtx_to_sce,
              mtxFile=snakemake@input$mtxs,
              bxFile=snakemake@input$bxs,
              txFile=snakemake@input$txs,
              sample_id=snakemake@params$sample_ids) %>%
    do.call(what=cbind)

################################################################################
                                        # Attach Annotations
################################################################################

colData(sce) %<>%
    as_tibble %>%
    left_join(df_cell_annots, by='cell_id') %>%
    set_rownames(.$cell_id) %>%
    DataFrame %>%
    `[`(colnames(sce),)

rowRanges(sce) <- gr_txs[rownames(sce), ]

if (has_row_annots) {
    rowData(sce) <- df_tx_annots[rownames(sce), ]
}

################################################################################
                                        # Export SCE
################################################################################

saveRDS(sce, snakemake@output$sce)
