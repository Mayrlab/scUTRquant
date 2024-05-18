#!/bin/bash

# mm10 UTRome 2022.08
wget -O tmp.tar.gz "https://figshare.com/ndownloader/files/36677958?private_link=78bbfa3c5e58a8cc5275" \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

# mm10 Annotations
## txs
wget -O utrome_mm10_v2_tx_annots.2024.05.13.csv.gz "https://figshare.com/ndownloader/files/46407730?private_link=78bbfa3c5e58a8cc5275"
wget -O utrome_mm10_v2_tx_annots.2024.05.13.Rds "https://figshare.com/ndownloader/files/46407733?private_link=78bbfa3c5e58a8cc5275"

## genes
wget -O utrome_mm10_v2_gene_annots.2024.05.13.csv.gz "https://figshare.com/ndownloader/files/46407724?private_link=78bbfa3c5e58a8cc5275"
wget -O utrome_mm10_v2_gene_annots.2024.05.13.Rds "https://figshare.com/ndownloader/files/46407727?private_link=78bbfa3c5e58a8cc5275"
