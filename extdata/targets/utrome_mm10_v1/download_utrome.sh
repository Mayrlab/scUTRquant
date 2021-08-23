#!/usr/bin/env bash

# mm10 UTRome 2019.05
wget -O tmp.tar.gz "https://ndownloader.figshare.com/files/24989696?private_link=2d0776b1d28138a5e716" \
    && tar -xvzf tmp.tar.gz --strip-components=1 \
    && rm tmp.tar.gz

# mm10 Atlas 2021.07
## txs
wget -O utrome_txs_annotation.Rds 'https://ndownloader.figshare.com/files/28873656?private_link=fa4677423a901fe137b2'

## genes
wget -O utrome_genes_annotation.Rds 'https://ndownloader.figshare.com/files/28873650?private_link=fa4677423a901fe137b2'
