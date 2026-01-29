#!/usr/bin/env bash


wget -O pbmc_unsorted_3k_gex.bam https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_gex_possorted_bam.bam

# Download analysis data and make annots.csv
wget -O tmp.tar.gz https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_analysis.tar.gz \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

join -t',' analysis/clustering/gex/graphclust/clusters.csv \
    analysis/dimensionality_reduction/gex/umap_projection.csv | \
    sed -E 's/([ACGT]{16})-1/pbmc_unsorted_3k_multiome_bam_\1/' | \
    sed '1 s/^.*$/cell_id,cluster,umap_1,umap_2/' > annots.csv
