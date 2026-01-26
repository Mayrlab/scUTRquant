#!/usr/bin/env bash

# mm10 UTRome 2019.05
wget -O tmp.tar.gz "https://github.com/Mayrlab/mca-utrome/releases/download/v1.0/utrome.w500.2019.05.tar.gz" \
    && tar -xvzf tmp.tar.gz --strip-components=1 \
    && rm tmp.tar.gz

# mm10 Atlas 2021.07
## txs
wget -O utrome_txs_annotation.Rds 'https://github.com/Mayrlab/atlas-mm/releases/download/v0.1/utrome_txs_annotation.Rds'

## genes
wget -O utrome_genes_annotation.Rds 'https://github.com/Mayrlab/atlas-mm/releases/download/v0.1/utrome_genes_annotation.Rds'
