class: middle, center

# Spring 2019 Reproducibility Workshop: "Containerization and Docker"

.left[Ryan Dale, Ph.D.]
.left[Scientific Information Officer, NICHD]
.left[Head, NICHD Bioinformatics and Scientific Programming Core]

Slides & code available at https://github.com/daler/docker-rw2019

---

class: middle, center

*Opinions are my own; nothing here should be construed as an endorsement*

---

# Containerization

Containerization is a way of packaging an application with all of its
dependencies into a standardized unit.

--

Helps resolve "well, it worked on *my* machine"

--

A container can be distributed and shared, ensuring others are using the same
programs

--

Popular containerization engines are Docker, rkt, Singularity.

---

# Docker

This document covers Docker, currently the most popular containerization
system.

A major limitation of Docker is that *it must run with root privileges*.

--

It is therefore not allowed on shared systems like HPC clusters. On NIH's
Biowulf cluster, use Singularity (https://hpc.nih.gov/apps/singularity.html).

---

# Docker, cloud, and scientific workflows

Web services use Docker as microservices on cloud instances listening on
different ports acting on incoming data.

--

This is a very different architecture from what is useful in science where we
have:

- many different interdependent tools
- interactive exploration
- scripts tied together
- heterogeneous data from many sources

--

Nevertheless, it is useful invest the effort to containerize especially at the
end of an analysis for long-term reproducibility.

--

**...but it's hard.**

---

# Terminology

**Docker Engine:** You need to install this on your computer to run Docker

--

**Container:** Standardized unit of software. Runs an OS, includes any
dependencies plus the application.

--

**Image:** A file which is a snapshot of a container. 

--

**Dockerfile:** Defines how to create an image. This is where you put your
effort when building a docker image

--

**Registry:** online location that stores Docker images. Docker Hub and quay.io
are the big ones.

--

**Tag:** Label for an image on a registry. By convention holds version information.

---

# Terminology usage

We write a *Dockerfile* which we use to build an *image*.

We run the image to create a running *container*.

Someone else can run that same image and get an identical running container.

---

# Get started

Prerequisites:

- [install Docker](https://docs.docker.com/v17.12/install/)
- Succesfully run the following command:

```bash
docker run hello-world
```

---

class: middle, center

# Worked example

---

# Overview of BioContainers

When working with Docker in bioinformatics, images for thousands of tools are
already available for your use.

This is a joint project between [bioconda](https://bioconda.github.io) and
[BioContainers](https://biocontainers.pro/#/).

We will use an existing container for `samtools`.

---

# Find an image to use


- Visit https://bioconda.github.io

- Search for `samtools`

- It tells us the image is `quay.io/biocontainers/samtools:<tag>` but we need to determine which tag to use.

- Follow the link on that site to visit https://quay.io/repository/biocontainers/samtools?tab=tags

- Use the latest tag unless you have a good reason for another version

--

So we want to use **`quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6`**

---

# Download an example BAM file

Download a file to use for demonstrating passing data into containers.

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

---

# Run the container

```bash
docker run quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6
```

--

**Explanation**

- download the image if needed
- run the image as a new container
- since there's no default command, after downloading and running nothing
  happened.

---

# Run the container again

We need to give it a command to run (`samtools`):

```bash
docker run quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 samtools
```

--

**Explanation**

- download the image if needed (it didn't need to this time)
- run the image as a new container
- within that running container, call `samtools`
- the help for `samtools` is printed to `stdout` (its default behavior
  when run with no arguments)

---

# Interactive containers

Use `-it` to drop into an interactive shell in the container.

**Anything you do here is ephemeral.**

```bash
docker run -it quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6
```

--

**Explanation:**

- run the image as a new container
- instead of exiting immediately, drop into an interactive shell
- type `exit` to exit the container

---

# View the downloaded file using the container

The BAM file we downloaded is in binary format; to read it we use `samtools`.

We might logically try the following.


```bash
# Note the "\" line continuation

docker run quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  samtools view y.bam
```

--

but this returns:

```
open: No such file or directory
[main_samview] fail to open "y.bam" for reading.
```

What happened?

---

# Make data available to the container

A running docker container is isolated. It cannot see anything on your
computer.

We need to explicitly tell it what is available with **`-v`**.

--

This allows us to mount paths from the host machine to paths inside the container.

```bash
docker run -v <host path>:<container path> <imagename>
```

--

**Notes:**

- Paths must be absolute
- Use **`$(pwd)`** as a shortcut for "current directory".
- Paths in the container are automatically created recursively
- `-v` should come *before* the image name, otherwise you'll get errors like
  this:

```text
docker: Error response from daemon: OCI runtime create failed:
container_linux.go:344: starting container process caused "exec: \"-v\":
executable file not found in $PATH": unknown.
```

---

# Mount the current directory to the container

Use `-v` to mount local paths inside the container:

```bash
docker run \
  -v $(pwd):/data quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  samtools view data/y.bam
```

--

Output should be:


```
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:1211:12450:44710     16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     HGIIFGEIJJJJJJJIJJJJIJIJJJJJIJJJJJJJ    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:1212:8027:29378      16      chr1    16286   255     3S33M   *       0       0    ATCTACACTGGGAGACACAGCAGTGAAGCTGAAATG     HG>DD<D?BD?@HGIGGGCFBCB<?9<IGIIIIGEF    NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0MD:Z:33 jM:B:c,-1       jI:B:i,-1       RG:Z:foo
AAGGGTCGT:DBR4KXP1:220:C27T3ACXX:7:1313:16438:83011     16      chr1    16286   255     3S33M   *       0       0    ATCTACAGTGGGAGACACAGCAGTGAAGCTGAAATG     B>GD?BB?)IEDBFFDHCEF?>GD9EEDGGCHFAGE    NH:i:1  HI:i:1  AS:i:30 nM:i:1  NM:i:1MD:Z:4C28       jM:B:c,-1       jI:B:i,-1       RG:Z:foo

...output truncated for brevity
```

---

# Mount the current directory to the container

Use `-v` to mount local paths inside the container:

```bash
docker run \
  -v $(pwd):/data quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  samtools view data/y.bam
```

**Explanation**

- `$(pwd)` is evaluated on the host system
- current working directory *outside* container becomes the directory `/data`
  *inside* the container
- `/data` was created automatically inside the container at run time
- The result is that the container can see everything in the working directory
  at `/data`, so the command run inside the container reflects this location
  and we get output as expected



---

# Get data back out of the container

`stdout` comes back to the host.

If the results you need are printed to `stdout`, you can capture results like
this:

```bash
docker run \
  -v $(pwd):/data \
  quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  samtools view /data/y.bam > result.sam
```

--

**Explanation**

- run `samtools` in the container, which prints to `stdout`
- bash on the host interprets **`>`** as the end of a command, and therefore the end of the `docker run` call
- `stdout` is redirected to a file *back on the host*, `result.sam`

---

# More complex input/output

Use `bash -c` to run more complex commands. This also demonstrates  mapping
multiple locations to the container:

```bash
mkdir -p results

docker run \
  -v $(pwd):/data \
  -v $(pwd)/results:/tmp \
  quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  bash -c "samtools view /data/y.bam > /tmp/result.sam"
```

--

**Explanation**

- after creating a local directory `results`, mount it to the container's `/tmp` directory
- `bash -c` means "commands are to be read from the following string"
- samtools output redirects to a file in `/tmp`, which is mapped to the local
  directory `results`.
- since the redirect `>` happens in the command executed by bash inside the
  container, we use container-appropriate paths.
- **user/group ownership of created files is `root`**

---

# Modify the user

```bash
mkdir -p results

docker run \
  -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  -v $(pwd)/results:/tmp \
  quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6 \
  bash -c "samtools view /data/y.bam > /tmp/result.sam"
```

**Explanation:**

- Run docker as the current user and group (*`id -u`* returns the UID of the
  current user)

---

# Combine multiple tools

This script downloads 100 reads from SRA and runs FastQC on the resulting file.

It uses a Docker container for sra-tools and another for FastQC.

This script is executed on the host, and communicates between containers by
mounting the same directory to each.

```bash
#!/bin/bash

# Download some FASTQ data using fastq-dump from sra-tools:
docker run \
  -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  quay.io/biocontainers/sra-tools:2.9.1_1--h470a237_0 \
  bash -c "cd /data && fastq-dump -X 100 --gzip SRR1574251"

# Run FastQC on the resulting file
docker run \
  -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  quay.io/biocontainers/fastqc:0.11.8--1 \
  bash -c "fastqc /data/*.fastq.gz"
```

---

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

---

# Build a custom container

BioContainers has containers for single tools.

What if we need multiple R packages in the same container?

We can't combine images; we need to build a custom image.

Use Bioconda to more easily build a custom image.

---

# Dockerfiles

A **Dockerfile** specifies how to build the container. The most
common Dockerfile commands are:

**`FROM`:** indicates which image should be used as a starting point. Here we use
`continuumio/miniconda3` because the conda package manager is already
installed. Generally `ubuntu:18.04` is a good choice.

**`RUN`:** Excutes the command in the building container.

**`COPY`:** Copy a file from the build directory on the host to a location
inside the image.

---

# `Dockerfile`

```docker
# This Dockerfile installs samtools, fastqc, R, and DESeq2

# We start from a container that already has conda installed
FROM continuumio/miniconda3

# Set up bioconda channel, see
# https://bioconda.github.io/index.html#set-up-channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Install packages. Pinning versions is optional but good practice
RUN conda install \
    "bioconductor-deseq2==1.22.1" \
    "fastqc==0.11.8" \
    "r-base>=3.5.1" \
    "samtools==1.9"
```

---

# Build the image

Tags provide easier-to-remember labels.

Conventional format is username/label.

Use your own Docker Hub username here if you're going to be uploading.

```bash
docker build . -t daler/rw2019
```

--

**Explanation**

- Execute the commands in the Dockerfile in the current directory (this may
  take some time)
- Save the resulting image as `daler/rw2019`

--

Optionally push to docker hub:

```bash
docker push daler/rw2019
```

---

# Put it all together

Run an R script, which is itself run by a bash script, all inside our new
container.

**Results:**

- FASTQ file
- FastQC report (HTML and zip)
- PCA plot in PDF format

---

# R script

Create an R script in the working directory.

It will be available in the container once the working directory is mounted.

```r
library(DESeq2)
dds <- makeExampleDESeqDataSet(betaSD=1)
rld <- rlog(dds)
pdf(file='pca.pdf')
plotPCA(rld)
dev.off()
sessionInfo()
```

---

# Modified bash script

Saved as `run2.sh`:

```bash
#!/bin/bash

cd /data
fastq-dump -X 100 --gzip SRR1574251
fastqc /data/*.fastq.gz
Rscript demo.R
```

---

# Run everything in a single container

```bash
docker run \
  -v $(pwd):/data \
  daler/rw2019 \
  /bin/bash run2.sh
```

**Explanation**

- in the container we built (which has everything we need), run the bash script
- the bash script also runs the R script

---

# Extra notes 1

- When building locally and then pushing, *Docker Hub does not show
  Dockerfiles* since the Dockerfile is not included in the image itself. It is
  better to have a GitHub repo and connect that to Docker Hub. Images will then
  get built automatically when the Dockerfile is changed on GitHub. Docker Hub
  will provide a link to that Dockerfile so that it can be inspected by the
  user.

- when building in Ubuntu containers, use the following idiom because apt repos
  are not included in the base images and therefore need to be updated. Be sure
  to include the `-y` argument, otherwise the build will hang waiting for user
  input.

```bash
RUN apt-get update && apt-get install -y \
  package1 \
  package2 \
```

---

# Extra notes 2

- Docker builds in layers, which are cached. Try to write Dockerfiles such that
  the earlier lines have less likelihood of needing to be changed.

- While it's possible to copy data over to the container, this wastes time and
  space -- try to maintain a separation of data and runtime dependencies

- Connect a GitHub repo to Docker Hub to auto-build containers (that users can
  see the Dockerfile for)

- One container -> one Snakemake rule

- One container -> entire environment for AWS instance
---

# Resources

- BioContainers [best-practices](https://biocontainers-edu.readthedocs.io/en/latest/best_practices.html)
- [BioContainers quickstart](http://biocontainers-edu.biocontainers.pro/en/latest/running_example.html)
- [Awesome Docker](https://github.com/veggiemonk/awesome-docker)
