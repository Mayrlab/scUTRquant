#!/usr/bin/env Rscript

library(magrittr)
library(tidyverse)

################################################################################
                                        # Fake Data (Interactive Testing)
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(
            ec=c("data/kallisto/utrome_mm10_v1/heart_1k_v2_fastq/matrix.ec"),
            txs=c("data/kallisto/utrome_mm10_v1/heart_1k_v2_fastq/transcripts.txt"),
            merge=c("data/utrs/utrome_mm10_v1/gene_merge.tsv")),
        output=list(tsv="/fscratch/fanslerm/test.ambiguity.tsv"),
        params=list())
}

################################################################################
                                        # Load Data
################################################################################

tx_map <- read_lines(snakemake@input$txs) %>%
    set_names(seq(0, along.with=.))

df_merge <- read_tsv(snakemake@input$merge,
                     col_names=c("tx_name", "merge_name"),
                     col_types="cc")

df_ec <- read_tsv(snakemake@input$ec,
                  col_names=c("ec", "tx"),
                  col_types="ic") %>%

    ## expand equivalence classes
    separate_rows(tx, sep=",", convert=FALSE) %>%

    ## map txs to tx_name
    mutate(tx_name=tx_map[tx]) %>%

    ## including merging table
    left_join(df_merge, by="tx_name") %>%

    ## does equivalence class include non-merging transcripts?
    group_by(ec) %>%
    mutate(is_multimapping=!all(first(merge_name) == merge_name)) %>%
    ungroup() %>%

    ## if so, call such txs "multimapping" and include num of ecs
    group_by(tx) %>%
    mutate(is_multimapping=any(is_multimapping), 
           n_ecs=n()) %>%

    ## does the transcript have ambiguity w/i the merging set?
    group_by(tx, merge_name) %>%
    mutate(is_overlapping=n() > 1) %>%
    ungroup() %>%

    ## retain tx-level summary
    select(merge_name, tx_name, n_ecs, is_multimapping, is_overlapping) %>%
    distinct()

################################################################################
                                        # Export table
################################################################################

write_tsv(df_ec, snakemake@output$tsv)
