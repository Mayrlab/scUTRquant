#!/usr/bin/env bash

# 10xv1
wget https://github.com/10XGenomics/cellranger/raw/3.0.2/lib/python/cellranger/barcodes/737K-april-2014_rc.txt

# 10xv2
wget https://github.com/10XGenomics/cellranger/raw/3.0.2/lib/python/cellranger/barcodes/737K-august-2016.txt

# 10xv3
wget https://github.com/10XGenomics/cellranger/raw/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt.gz \
    && gzip -d 3M-february-2018.txt.gz

