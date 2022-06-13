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
            bxs=c("data/kallisto/neuron_1k_v2_fastq/utrome.genes.barcodes.txt"),
            genes=c("data/kallisto/neuron_1k_v2_fastq/utrome.genes.genes.txt"),
            mtxs=c("data/kallisto/neuron_1k_v2_fastq/utrome.genes.mtx"),
            gtf="extdata/gff/adult.utrome.e3.t200.f0.999.w500.gtf",
            gene_annots="extdata/atlas/utrome_genes_annotation.Rds"),
        output=list(sce="data/sce/test.utrome.genes.Rds"),
        params=list(genome='mm10',
                    sample_ids=c("neuron_1k_v2_fastq"),
                    min_umis="500",
                    cell_annots="examples/neuron_1k_v2_fastq/annots.csv"))
}

################################################################################
                                        # Load Data
################################################################################

## Row Annotations (UTRs)
has_row_annots <- !is.null(snakemake@input$gene_annots)
if (has_row_annots) {
    df_gene_annots <- readRDS(snakemake@input$gene_annots)
}

## Row Ranges (
gr_genes <- read_gff(snakemake@input$gtf, genome_info=snakemake@params$genome,
                    col_names=c('type', 'transcript_id', 'gene_id')) %>%

    ## keep only transcripts
    filter(type == 'transcript') %>%
    select(-c('type')) %>%

    ## add rownames
    `names<-`(.$transcript_id) %>%

    ## group by gene_id
    split(.$gene_id)

## Column Annotations
df_cell_annots <- snakemake@params$cell_annots %>% {
    if (is.null(.)) {
        tibble(cell_id=character(0))
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
                               colData=DataFrame(cell_id=colnames(.), sample_id=sample_id, row.names=colnames(.))) }
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
    left_join(df_cell_annots, by='cell_id') %>%
    set_rownames(.$cell_id) %>%
    DataFrame %>%
    `[`(colnames(sce),)


## TODO: Change merge table to use gene_id upstream
## switch to gene_id instead of gene_symbol
if (has_row_annots) {
    sym2id <- df_gene_annots[,c("gene_symbol", "gene_id")] %>% deframe()
    rownames(sce) <- sym2id[rownames(sce)]
}

## attach location info
rowRanges(sce) <- gr_genes[rownames(sce), ]

if (has_row_annots) {
    rowData(sce) <- df_gene_annots[rownames(sce), ]
}

################################################################################
                                        # Export SCE
################################################################################

saveRDS(sce, snakemake@output$sce)
