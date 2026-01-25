#!/bin/bash

# mm10 UTRome 2022.08
wget -O tmp.tar.gz "https://github.com/Mayrlab/mca-utrome/releases/download/v2.0.1/utrome.e30.t5.gc25.pas3.f0.9999.w500.tar.gz" \
    && tar -xvzf tmp.tar.gz \
    && rm tmp.tar.gz

# mm10 Annotations
## txs
wget -O utrome_mm10_v2_tx_annots.2024.05.13.csv.gz "https://github.com/Mayrlab/atlas-mm/releases/download/v0.2.0/utrome_mm10_v2_tx_annots.2024.05.13.csv.gz"
wget -O utrome_mm10_v2_tx_annots.2024.05.13.Rds "https://github.com/Mayrlab/atlas-mm/releases/download/v0.2.0/utrome_mm10_v2_tx_annots.2024.05.13.Rds"

## genes
wget -O utrome_mm10_v2_gene_annots.2024.05.13.csv.gz "https://github.com/Mayrlab/atlas-mm/releases/download/v0.2.0/utrome_mm10_v2_gene_annots.2024.05.13.csv.gz"
wget -O utrome_mm10_v2_gene_annots.2024.05.13.Rds "https://github.com/Mayrlab/atlas-mm/releases/download/v0.2.0/utrome_mm10_v2_gene_annots.2024.05.13.Rds"
