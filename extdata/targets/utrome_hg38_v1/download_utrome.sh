#!/bin/bash

# hg38 UTRome 2022.06
wget -O tmp.tar.gz "https://figshare.com/ndownloader/files/35874803?private_link=9f97fd9772574bee6a61" \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

# hg38 Annotations
## txs
wget -O utrome_hg38_v1_tx_annots.2024.05.13.csv.gz "https://figshare.com/ndownloader/files/46407826?private_link=9f97fd9772574bee6a61"
wget -O utrome_hg38_v1_tx_annots.2024.05.13.Rds "https://figshare.com/ndownloader/files/46407829?private_link=9f97fd9772574bee6a61"

## genes
wget -O utrome_hg38_v1_gene_annots.2024.05.13.csv.gz "https://figshare.com/ndownloader/files/46407820?private_link=9f97fd9772574bee6a61"
wget -O utrome_hg38_v1_gene_annots.2024.05.13.Rds "https://figshare.com/ndownloader/files/46407823?private_link=9f97fd9772574bee6a61"
