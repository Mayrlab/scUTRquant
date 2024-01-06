#!/usr/bin/env bash

# 10xv1
wget https://teichlab.github.io/scg_lib_structs/data/737K-april-2014_rc.txt.gz \
     && gzip -d 737K-april-2014_rc.txt.gz

# 10xv2
wget https://teichlab.github.io/scg_lib_structs/data/737K-august-2016.txt.gz \
     && gzip -d 737K-august-2016.txt.gz

# 10xv3
wget https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz \
    && gzip -d 3M-february-2018.txt.gz
