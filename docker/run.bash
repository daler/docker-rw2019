#!/bin/bash
cd /data
fastq-dump -X 100 --gzip SRR1574251
fastqc /data/*.fastq.gz
Rscript demo.R > R.log
