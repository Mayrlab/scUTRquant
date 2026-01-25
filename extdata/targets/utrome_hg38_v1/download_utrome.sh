#!/bin/bash

# hg38 UTRome 2022.06
wget -O tmp.tar.gz "https://github.com/Mayrlab/hcl-utrome/releases/download/v1.0.0/utrome.e30.t5.gc39.pas3.f0.9999.w500.tar.gz" \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

# hg38 Annotations
## txs
wget -O utrome_hg38_v1_tx_annots.2024.05.13.csv.gz "https://github.com/Mayrlab/atlas-hs/releases/download/v0.1.0/utrome_hg38_v1_tx_annots.2024.05.13.csv.gz"
wget -O utrome_hg38_v1_tx_annots.2024.05.13.Rds "https://github.com/Mayrlab/atlas-hs/releases/download/v0.1.0/utrome_hg38_v1_tx_annots.2024.05.13.Rds"

## genes
wget -O utrome_hg38_v1_gene_annots.2024.05.13.csv.gz "https://github.com/Mayrlab/atlas-hs/releases/download/v0.1.0/utrome_hg38_v1_gene_annots.2024.05.13.csv.gz"
wget -O utrome_hg38_v1_gene_annots.2024.05.13.Rds "https://github.com/Mayrlab/atlas-hs/releases/download/v0.1.0/utrome_hg38_v1_gene_annots.2024.05.13.Rds"
