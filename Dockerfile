FROM mambaorg/micromamba:latest

LABEL maintainer="matt@datatecnica.com"
LABEL version="1.0.0"
LABEL image="colocalization"
LABEL description="Container image for QTL + GWAS colocalization analysis pipeline."

USER root

# Install system dependencies required for R packages, bioinformatics tools, and compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    wget \
    libgl1 \
    libglib2.0-0 \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libudunits2-dev \
    libgdal-dev \
    libmagick++-dev \
    librsvg2-dev \
    libfftw3-dev \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    libdeflate-dev \
    libuv1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy environment definition and create the conda environment
COPY env/environment.yml /tmp/environment.yml
RUN micromamba env create -f /tmp/environment.yml \
    && micromamba clean -afy

# Activate the colocalization environment for all subsequent RUN steps
SHELL ["micromamba", "run", "-n", "colocalization", "/bin/bash", "-c"]

# Copy and run the R library installation script
COPY env/install_libraries.R /tmp/install_libraries.R
RUN Rscript /tmp/install_libraries.R

# Copy the full project into the container
WORKDIR /colocalization
COPY . .

# Ensure the colocalization environment is on PATH at runtime
ENV PATH="/opt/conda/envs/colocalization/bin:$PATH"
ENV CONDA_DEFAULT_ENV=colocalization

# Default: run the GTEx analysis workflow; override ENTRYPOINT at runtime for individual steps
ENTRYPOINT ["micromamba", "run", "-n", "colocalization", "/bin/bash"]
