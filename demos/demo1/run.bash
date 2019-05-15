#!/bin/bash

# Download some FASTQ data using fastq-dump from sra-tools:
docker run \
  -v $(pwd):/data \
  quay.io/biocontainers/sra-tools:2.9.1_1--h470a237_0 \
  bash -c "cd /data && fastq-dump -X 100 --gzip SRR1574251"

# Run FastQC on the resulting file
docker run \
  -v $(pwd):/data \
  quay.io/biocontainers/fastqc:0.11.8--1 \
  bash -c "fastqc /data/*.fastq.gz"
