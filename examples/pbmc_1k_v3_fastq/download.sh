#!/usr/bin/env bash

wget -O tmp.tar https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar \
  && tar -xf tmp.tar \
  && mv pbmc_1k_v3_fastqs fastq \
  && rm tmp.tar

## Download analysis data and make annots.csv
wget -O tmp.tar.gz https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_analysis.tar.gz \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

join -t',' analysis/clustering/graphclust/clusters.csv analysis/tsne/2_components/projection.csv | \
    sed -E 's/([ACGT]{16})-1/pbmc_1k_v3_fastq_\1/' | \
    sed '1 s/^.*$/cell_id,cluster,tsne_1,tsne_2/' > annots.csv
