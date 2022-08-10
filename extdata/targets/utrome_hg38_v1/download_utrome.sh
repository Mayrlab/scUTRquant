#!/bin/bash

# hg38 UTRome 2022.06
wget -O tmp.tar.gz "https://figshare.com/ndownloader/files/35874803?private_link=9f97fd9772574bee6a61" \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

# hg38 Annotations
## txs
#wget -O utrome_txs_annotation.Rds

## genes
#wget -O utrome_genes_annotation.Rds
