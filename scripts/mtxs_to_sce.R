library(magrittr)
library(plyranges)
library(tidyverse)
library(Matrix)
library(SingleCellExperiment)

################################################################################
                                        # Fake Data (Interactive Testing)
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(
            bxs=c("data/kallisto/utrome/lsk_sigab1/utrome.txs.barcodes.txt",
                  "data/kallisto/utrome/lk_sigaf1/utrome.txs.barcodes.txt"),
            txs=c("data/kallisto/utrome/lsk_sigab1/utrome.txs.genes.txt",
                  "data/kallisto/utrome/lk_sigaf1/utrome.txs.genes.txt"),
            mtxs=c("data/kallisto/utrome/lsk_sigab1/utrome.txs.mtx",
                   "data/kallisto/utrome/lk_sigaf1/utrome.txs.mtx"),
            gtf="/data/mayrc/data/mca/gff/adult.utrome.e3.t200.f0.999.w500.gtf",
            annots="metadata/dahlin18_wolf19_annots.csv"),
        output=list(sce="/fscratch/fanslerm/hspcs.utrome.txs.Rds"),
        params=list(genome='mm10',
                    sample_ids=c("SIGAB1", "SIGAF1"),
                    min_umis="1000"))
}

################################################################################
                                        # Load Data
################################################################################

## Row Annotations (UTRs)
gr_utrs <- read_gff(snakemake@input$gtf, genome_info=snakemake@params$genome,
                    col_names=c('type', 'transcript_id', 'transcript_name', 'gene_id')) %>%

    ## keep only transcripts
    filter(type == 'transcript') %>%
    select(-c('type')) %>%

    ## add gene_symbol
    mutate(gene_symbol=str_remove(transcript_id, "\\.[^.]+$")) %>%

    ## add rownames
    `names<-`(.$transcript_id)

## Column Annotations
## df_annots <- read_csv(snakemake@input$annots)

## Load MTXs Function
load_mtx_to_sce <- function (mtxFile, bxFile, txFile, sample_id) {
    ## TODO: Add annotation-based filtering
    ## valid_cells <- filter(df_annots, sample_id == sample_id) %$% cell_id
    readMM(mtxFile) %>%
        as("CsparseMatrix") %>%
        `dimnames<-`(list(
            cell_id=str_c(read_lines(bxFile), "_", sample_id),
            transcript_id=read_lines(txFile))) %>%
        t %>%
        { .[, colSums(.) >= as.integer(snakemake@params$min_umis)] } %>%
        { SingleCellExperiment(assays=list(counts=.)) }
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

## colData(sce) <- df_annots %>%
##     set_rownames(.$cell_id) %>%
##     DataFrame %>%
##     `[`(colnames(sce),)

rowRanges(sce) <- gr_utrs[rownames(sce), ]

################################################################################
                                        # Export SCE
################################################################################

saveRDS(sce, snakemake@output$sce)
