#!/usr/bin/env #!/usr/bin/env bash

wget -O tmp.tar https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_1k_v2/heart_1k_v2_fastqs.tar \
  && tar -xf tmp.tar \
  && mv heart_1k_v2_fastqs fastq \
  && rm tmp.tar
