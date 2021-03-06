# hybrid-assembly
# UNDER DEVELOPMENT

[![Build Status](https://travis-ci.org/kevinmenden/hybrid-assembly.svg?branch=master)](https://travis-ci.org/kevinmenden/hybrid-assembly)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Repository on Dockerhub](https://img.shields.io/badge/docker-available-green.svg "Docker Repository on Dockerhub")](https://hub.docker.com/r/kevinmenden/hybrid-assembly/)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The hybrid-assembly pipeline contains workflows for de-novo genome assembly using both short and long reads.
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

The current workflow consists of:
1. Quality control of input reads
2. Assembly
3. Assesment of assembly quality

Two different assemblers are currently implemented:
* [FLYE](https://github.com/fenderglass/Flye)
* [MaSuRCA](https://github.com/alekseyzimin/masurca)


### Documentation
The hybrid-assembly pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
This pipeline was written by Kevin Menden ([kevinmenden](https://github.com/kevinmenden)) at [DZNE](http://www.dzne.de).

[![DZNE](assets/dzne-logo.jpeg)](http://www.dzne.de)
