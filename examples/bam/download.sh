#!/usr/bin/env bash

wget -O neuron_1k_v3.bam https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_possorted_genome_bam.bam

wget -qO- https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_analysis.tar.gz | tar xvz -
