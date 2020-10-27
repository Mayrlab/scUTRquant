#!/usr/bin/env bash

wget -O neuron_1k_v3.bam https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_possorted_genome_bam.bam

# Download analysis data and make annots.csv
# wget -O tmp.tar.gz https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_analysis.tar.gz \
#     && tar -xvzf tmp.tar.gz \
#     && rm tmp.tar.gz

# join -t',' analysis/clustering/graphclust/clusters.csv analysis/tsne/2_components/projection.csv | \
#     sed -E 's/([ACGT]{16})-1/neuron_1k_v3_bam_\1/' | \
#     sed '1 s/^.*$/cell_id,cluster,tsne_1,tsne_2/' > annots.csv
