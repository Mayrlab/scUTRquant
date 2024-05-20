FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8

# Install required packages in base
COPY --chown=$MAMBA_USER:$MAMBA_USER scutr-quant.yaml /tmp/scutr-quant.yaml
RUN micromamba install -yn base -f /tmp/scutr-quant.yaml \
  && micromamba clean -qafy

# Put base bin on PATH
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
