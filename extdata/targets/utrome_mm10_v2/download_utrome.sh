#!/bin/bash

# mm10 UTRome 2022.08
wget -O tmp.tar.gz "https://figshare.com/ndownloader/files/36677958?private_link=78bbfa3c5e58a8cc5275" \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

# mm10 Annotations
## txs
wget -O utrome_txs_annotation.Rds "https://figshare.com/ndownloader/files/36741822?private_link=78bbfa3c5e58a8cc5275"

## genes
wget -O utrome_genes_annotation.Rds "https://figshare.com/ndownloader/files/36741819?private_link=78bbfa3c5e58a8cc5275"
