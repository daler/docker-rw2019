# Spring 2019 Reproducibility Workshop: "Containerization and Docker"

Ryan Dale, Ph.D.

*Opinions are my own; nothing here should be construed as an endorsement*

## Overview

Containerization is a way of packaging an application with all of its
dependencies into a standardized unit.

This helps solve the "well it works on *my* machine" problem of
reproducibility.

Popular containerization engines are Docker, rkt, Singularity.

This document covers Docker. A severe limitation of Docker is that *it must run
with root privileges*. It is therefore not allowed on shared systems like HPC
clusters.

On NIH's Biowulf cluster, use Singularity
(https://hpc.nih.gov/apps/singularity.html). It plays nicely with Docker
images.

## Terminology

**Docker Engine:** You need to install this on your computer to run containers
or build images

**Container:** standardized unit of software. Runs an OS, includes any
dependencies plus the application

**Image:** holds layers that, when mounted on the Docker engine, operate as
a container

**Dockerfile:** Defines how to create an image. This is where you put your
effort when building a docker image. Can inherit from other images, install
software, run things by default

**Registry:** online location that stores Docker images. Docker Hub and quay.io
are the big ones.

**Tag:** Label for an image on a registry. By convention holds version information.

## Getting started

Prerequisites:

- [install Docker](https://docs.docker.com/v17.12/install/)
- Succesfully run the following command:

```bash
docker run -it hello-world
```

## Worked example

This example gradually builds up in complexity.

### Download an example SAM file

Let's get some example data to work with. This is a small SAM file from the
samtools GitHub repo that we'll use to demonstrate how to use Docker.

With `wget`:

```bash
wget https://github.com/daler/pybedtools/raw/master/pybedtools/test/data/y.bam
```

or with `curl`:

```bash
curl https://github.com/daler/pybedtools/raw/master/pybedtools/test/data/y.bam > y.bam
```

Or download the file with web browser from:
https://github.com/daler/pybedtools/raw/master/pybedtools/test/data/y.bam


### Overview of BioContainers

When working with Docker in bioinformatics, thousands of images are already
available for your use. We will use an existing container for `samtools`. Here
is the process of finding which image to use:

- Visit https://bioconda.github.io
- Search for `samtools`
- It tells us the image is `quay.io/biocontainers/samtools:<tag>` but we need
  to visit https://quay.io/repository/biocontainers/samtools?tab=tags to see
  our options are for tags
- Use the latest tag unless you have a good reason for another version

So we want to use `quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6`.

Here's how to use it. This will pull and run the container:

```bash
docker run quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6
```

Since there's no default command, nothing happened. So we specify the command
to run (`samtools`):

```bash
docker run quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 samtools
```

Explanation:

- download the image if needed
- run the image as a new container
- within that running container, call `samtools`
- the help for `samtools` is printed to `stdout` which is its default behavior
  when run with no arguments

### Try viewing the downloaded file

To view the reads in that SAM file we downloaded using `samtools`, we might
logically try the following:

```bash
docker run quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 samtools view y.bam
```

but this returns:

```
[main_samview] fail to open "y.bam" for reading.
```

What happened?

### Making data available to the container

A running docker container is isolated. It cannot see anything on your
computer; you need to explicitly tell it what is available.

Use `-v` to mount paths from your local machine to paths inside the container.

- Paths must be absolute, so we typically use `$(pwd)` as a shortcut for "current directory".
- Paths in the container will be created recursively as needed
- The `-v` command should come before the image name, otherwise you'll get errors like this:

```
docker: Error response from daemon: OCI runtime create failed:
container_linux.go:344: starting container process caused "exec: \"-v\":
executable file not found in $PATH": unknown.
```

This works:

```bash
docker run \
  -v $(pwd):/data quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  samtools view data/y.bam
```

Output should be:

```
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:1211:12450:44710     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     HGIIFGEIJJJJJJJIJJJJIJIJJJJJIJJJJJJJ    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:1212:8027:29378      16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     HG>DD<D?BD?@HGIGGGCFBCB<?9<IGIIIIGEF    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:1313:16438:83011     16      chr1    16286   255     3S33M   *       0       0    ATCTACAGTGGGAGACACAGCAGTGAAGCTGAAATG     B>GD?BB?)IEDBFFDHCEF?>GD9EEDGGCHFAGE    NH:i:1  HI:i:1  AS:i:30 nM:i:1  NM:i:1MD:Z:4C28       jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2105:9469:36427      16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     GD<?<AFGEIHGFCGCGGFFDGGAEGFFBEHGHHDG    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2108:17666:14770     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     HCEIGIHJJJJIGIHGJIIJJJJJJJJIJJIJJJIJ    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2109:17746:60505     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     IDHFCGJJJIIIIIGJJJIJJJIGJJIHDIJHEGIH    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2109:4668:57747      16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     IDIIHEHGBB@IJJIIHGGIGGIGGJGIGHBJIHHD    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2207:15896:13562     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     IHHEIIHJJJJJJIJJJIIIHGJJIIGHEIJJIIGI    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2210:15325:55969     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     HCDDB9?GIGFBBDHGGEGIHHHCGEABGCCDDC4+    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2211:18156:35604     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     IGFD?@FADEGGEIIGF@CEGEFGDHGGEJIJIJIJ    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2301:14313:83259     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     HD9@GDDGDHGDEFFCGGDB<CF>DHGFE@AIHGC?    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:2306:1697:8254       16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     JJJJIJJJIJJJJJJJJJJJJJJJJJJJIJIJJJJJ    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:1302:9624:39969      16      chr1    16292   255     27M     *       0       0    GGGAGACACAGCAGTGAAGCTGAAATG      DB???:0*:1**9:*?C4?119FC9<>     NH:i:1  HI:i:1  AS:i:26 nM:i:0  NM:i:0  MD:Z:27 jM:B:c,-1     jI:B:i,-1       RG:Z:foo
TTCGGTCAC:DBR4KXP1:220:C27T3ACXX:7:1109:13559:30319     16      chr1    16439   255     41M     *       0       0    CTCTACAGTTTGAAAACCACTATTTTATGAACCAAGTAGAA        @GBB<CBDF@??GEDGDC:<FF@IHBHD@@AA<EF>HAFAH       NH:i:1  HI:i:1  AS:i:40       nM:i:0  NM:i:0  MD:Z:41 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
TTCGGTCAC:DBR4KXP1:220:C27T3ACXX:7:1305:1518:19706      16      chr1    16439   255     41M     *       0       0    CTCTACAGTTTGAAAACCACTATTTTATGAACCAAGTAGAA        @GIIIGGIIGGIGFE@GFEIGGIHCGIEIHHIIIIHGBFFH       NH:i:1  HI:i:1  AS:i:40       nM:i:0  NM:i:0  MD:Z:41 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
TTCGGTCAC:DBR4KXP1:220:C27T3ACXX:7:1306:17189:41340     16      chr1    16439   255     41M     *       0       0    CTCTACAGTTTGAAAACCACTATTTTATGAACCAAGTAGAA        GJIIGCEGIIIGIFEGEFBHEIJIIGIGIGHGHGGJIDHFH       NH:i:1  HI:i:1  AS:i:40       nM:i:0  NM:i:0  MD:Z:41 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
TTCGGTCAC:DBR4KXP1:220:C27T3ACXX:7:2105:17701:25727     16      chr1    16439   255     41M     *       0       0    CTCTACAGTTTGAAAACCACTATTTTATGAACCAAGTAGAA        IJIIFIJJIJJJJJIGGE?JJJJJJJJJJIHIJJJJJHHHH       NH:i:1  HI:i:1  AS:i:40       nM:i:0  NM:i:0  MD:Z:41 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
TTCGGTCAC:DBR4KXP1:220:C27T3ACXX:7:2106:5853:94674      16      chr1    16439   255     41M     *       0       0    CTCTACAGTTTGAAAACCACTATTTTATGAACCAAGTAGAA        F?*004GD<GFFFBECC:1A>@DHFAIHF<:EEGFA4:<8D       NH:i:1  HI:i:1  AS:i:40       nM:i:0  NM:i:0  MD:Z:41 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
GCAGGTCAA:DBR4KXP1:220:C27T3ACXX:7:2113:10820:21294     0       chr1    30891   255     40M1S   *       0       0    TCATTTCTCTCTATCTCATTTCTCTCTCTCTCGCTCTTTCT        BBFAFFHIIGBH:CEFHIIIIIIIIHIHIIII)?GGHEHHI       NH:i:1  HI:i:1  AS:i:35       nM:i:2  NM:i:2  MD:Z:35A1C2     jM:B:c,-1       jI:B:i,-1       RG:Z:foo
GCCGGTCTC:DBR4KXP1:220:C27T3ACXX:7:1105:10358:96045     16      chr1    46637   255     41M     *       0       0    TCTGCCTTGAAATTCTTAACAATTTTTTTAACCAAAGTCCT        ;@==BFGGEFIIFHBGGHH@GGGGGGDGHFAG>D@FGFFFA       NH:i:1  HI:i:1  AS:i:40       nM:i:0  NM:i:0  MD:Z:41 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
GCCGGTCTC:DBR4KXP1:220:C27T3ACXX:7:1107:16157:39614     16      chr1    46637   255     41M     *       0       0    TCTGCCTTGAAATTCTTAACAATTTTTTTAACCAAAGTCCT        F;8;<B?@HCIHGB?FGD:<HEBE>HEC<?<@GIFGHHDFB       NH:i:1  HI:i:1  AS:i:40       nM:i:0  NM:i:0  MD:Z:41 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
```


Explanation:

- mount the current working directory *outside* the container to the directory
  `/data` *inside* the container
- `$(pwd)` is an argument evaluated on the host system, not in the docker
  container
- `/data` was created automatically inside the container at run time
- The result is that the container can see everything in the working directory
  at `/data`, so the command run inside the container reflects this location
  and we get output as expected

### Getting data back out of the container

`stdout` comes back to the host. So if the results you need are printed to
`stdout`, you can capture results like this:

```bash
docker run \
  -v $(pwd):/data \
  quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  samtools view /data/y.bam \
  > result.sam
```

Explanation:

- redirect the results to `result.sam` *back on the host*, not in the Docker
  container

### More complex input/output

This expands the last example to show how we can run more complex commands
using `bash -c` as well as mapping multiple locations.

```bash
mkdir -p results

docker run \
  -v $(pwd):/data \
  -v $(pwd)/results:/tmp \
  quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  bash -c "samtools view /data/y.bam > /tmp/result.sam
```

Explanation:

- after creating a local directory `results`, mount it to the container's `/tmp` directory
- `bash -c` means "commands are to be read from the following string"
- samtools output redirects to a file in `/tmp`, which is mapped to the local
  directory `results`.
- since the redirect `>` happens in the command executed by bash inside the
  container, we use container-appropriate paths.


## Combining multiple tools

Expanding the example further, this demonstrates how to run multiple tools.
This downloads 100 reads from SRA and runs FastQC on the resulting file.

This is a bash script that is executed on the host. It runs two separate
containers. Even though the containers are isolated from each other and the
host, by providing consistent `-v` argument we get them to see the same data.

```bash
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
```

Assuming the above file is saved as `run.sh`, run this locally as:

```bash
bash run.sh
```

You should get the following output:

```
Read 100 spots for SRR1574251
Written 100 spots for SRR1574251
Started analysis of SRR1574251.fastq.gz
Analysis complete for SRR1574251.fastq.gz
```

And when you check your working directory, you should see:

```
SRR1574251.fastq.gz
SRR1574251_fastqc.zip
SRR1574251_fastqc.html
```

## Building a custom container

This next example creates a larger container with R, DESeq2, FastQC, and
samtools. We use bioconda to save LOTS of time. This has the side effect of
being able to specify specific versions into the container.

We write a Dockerfile which specifies how to build the container. The most
common Dockerfile commands are **`FROM`** and **`RUN`**

**FROM:** indicates which image should be used as a starting point. Here we use
`continuumio/miniconda3` because the conda package manager is already
installed. Generally `ubuntu:18.04` is a good choice.

**RUN:** Excutes the command in the building container.

Some notes:

- when building in Ubuntu containers, use the following idiom because apt repos
  are not included in the base images and therefore need to be updated:

```bash
RUN apt-get update && apt-get install \
  package1 \
  package2 \
```

Here's our Dockerfile:

```
# This Dockerfile installs samtools, fastqc, R, and DESeq2

# We start from a container that already has conda installed
FROM continuumio/miniconda3

# Set up bioconda channel, see
# https://bioconda.github.io/index.html#set-up-channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Install packages
RUN conda install \
    "bioconductor-deseq2==1.22.1" \
    "fastqc==0.11.8" \
    "r-base>=3.5.1" \
    "samtools==1.9"
```

```bash
# Format is username/label. Use your own Docker Hub username here if you're
# going to be uploading.
docker build . -t daler/rw2019
```

Push to docker hub:

```bash
docker push daler/rw2019
```

Create an R script in the working directory. When mountin the current
directory, it will be accessible inside the container:

```r
library(DESeq2)
dds <- makeExampleDESeqDataSet(betaSD=1)
rld <- rlog(dds)
pdf(file='pca.pdf')
plotPCA(rld)
dev.off()
sessionInfo()
```

```bash
docker run \
  -v $(pwd):/data \
  daler/demo-container \
  Rscript demo.R
```

## Extra notes

When building locally and then pushing, *Docker Hub does not show Dockerfiles*
since the Dockerfile is not included in the image itself. It is better to have
a GitHub repo and connect that to Docker Hub. Images will then get built
automatically when the Dockerfile is changed on GitHub. Docker Hub will provide
a link to that Dockerfile so that it can be inspected by the user.


## Resources

- BioContainers [best-practices](https://biocontainers-edu.readthedocs.io/en/latest/best_practices.html)
- [BioContainers quickstart](http://biocontainers-edu.biocontainers.pro/en/latest/running_example.html)
- [Awesome Docker](https://github.com/veggiemonk/awesome-docker)
