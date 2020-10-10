#!/usr/bin/env #!/usr/bin/env bash

wget -O tmp.tar https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_fastqs.tar \
  && tar -xvzf tmp.tar \
  && mv neuron_1k_v3_fastqs fastq \
  && rm tmp.tar
