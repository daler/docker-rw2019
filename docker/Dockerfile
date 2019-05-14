# This Dockerfile installs samtools, fastqc, R, and DESeq2
# We start from a container that already has conda installed
FROM continuumio/miniconda3

# Set up bioconda channel, see
# https://bioconda.github.io/index.html#set-up-channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Install packages. Pinning versions is optional but good practice
#
RUN conda install \
    "bioconductor-deseq2==1.22.1" \
    "fastqc==0.11.8" \
    "r-base>=3.5.1" \
    "r-ggplot2" \
    "samtools==1.9"

# Alternative method, if requirements are stored outside the Dockerfile:
#
# COPY requirements.txt .
# RUN conda install --file requirements.txt
