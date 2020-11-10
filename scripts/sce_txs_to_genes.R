#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(Matrix)
library(magrittr)

################################################################################
                                        # Fake Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(sce="data/sce/neuron_10k_v3_fastq.utrome.txs.Rds"),
        output=list(sce="/fscratch/fanslerm/neuron_10k_v3_fastq.utrome.genes.Rds"),
        params=list())
}

################################################################################
                                        # Load SCE
################################################################################

sce <- readRDS(snakemake@input$sce)

################################################################################
                                        # Gene Counts
################################################################################

cts_genes <- fac2sparse(rowData(sce)$gene_id) %*% counts(sce)

gr_genes <- rowRanges(sce) %>% split(.$gene_id)

sce_genes <- SingleCellExperiment(assays=list(counts=cts_genes),
                                  rowRanges=gr_genes[rownames(cts_genes)],
                                  colData=colData(sce))

################################################################################
                                        # Export Data
################################################################################

saveRDS(sce_genes, snakemake@output$sce)
