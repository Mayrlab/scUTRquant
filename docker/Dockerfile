FROM continuumio/miniconda3

SHELL ["/bin/bash", "--login", "-c"]

# Install required packages in base
COPY scutr-quant.yaml .
RUN conda env update -n base -f scutr-quant.yaml \
  && conda clean -a -y

# Install custom kallisto build
COPY kallisto-build.yaml .
RUN conda env create -f kallisto-build.yaml -n kallisto-build \
  && git clone -b stranded --single-branch https://github.com/mfansler/kallisto.git /tmp/kallisto \
  && conda activate kallisto-build \
  && mkdir /tmp/kallisto/build && cd /tmp/kallisto/build \
  && conda activate kallisto-build \
  && cmake /tmp/kallisto && make && make install \
  && cd $HOME && rm -rf /tmp/kallisto \
  && conda deactivate \
  && conda env remove -n kallisto-build -y \
  && conda clean -a -y
