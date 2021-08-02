FROM --platform=linux/amd64 condaforge/mambaforge:4.10.3-3

SHELL ["/bin/bash", "--login", "-c"]

# Install custom kallisto build
COPY kallisto-build.yaml .
RUN mamba env create -f kallisto-build.yaml -n kallisto-build \
  && git clone -b stranded --single-branch https://github.com/mfansler/kallisto.git /tmp/kallisto \
  && mkdir -p /tmp/kallisto/build && cd /tmp/kallisto/build \
  && . /opt/conda/etc/profile.d/conda.sh \
  && conda activate kallisto-build \
  && cmake /tmp/kallisto && make && make install \
  && cd $HOME && rm -rf /tmp/kallisto \
  && conda deactivate \
  && mamba env remove -n kallisto-build -y \
  && mamba clean -qafy

# readr needs timezone info
RUN apt-get update > /dev/null \
  && DEBIAN_FRONTEND="noninteractive" apt-get install --no-install-recommends -y tzdata > /dev/null \
  && apt-get clean

# Install required packages in base
COPY scutr-quant.yaml .
RUN mamba env update -n base -f scutr-quant.yaml \
  && mamba clean -qafy
