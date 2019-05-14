# Spring 2019 Reproducibility Workshop

Download
[remark.js slides](https://raw.githubusercontent.com/daler/docker-rw2019/master/slides/slides.html)
and view locally in a browser.

The slides are created from a Markdown version 
([slides/slides.md](slides/slides.md)) using [slides/prep.py](slides/prep.py).
It may be easier to copy/paste from the Markdown file.

The [docker](docker) directory contains a Dockerfile and other configuration options.

Docker Hub has been connected to this GitHub repository, so that any changes to
[`docker/Dockerfile`](docker/Dockerfile) will result in a new image.

You can use the latest image like this:

```bash
docker run -it daler/docker-rw2019:latest
```

You can view the Docker Hub page for this repo at
https://hub.docker.com/r/daler/docker-rw2019.
