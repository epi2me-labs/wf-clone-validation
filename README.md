# Clone validation workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
that can be used to de novo assemble plasmid sequences from Oxford Nanopore data.

## Introduction

This repository contains a nextflow workflow that can be used to de novo assemble plasmid sequences from Oxford Nanopore data.
## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-clone-validation --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a consensus fasta file for each sample
* a csv showing the pass or fail status of each sample
* a feature table containing annotations for each of the samples
* an HTML report document detailing the primary findings of the workflow.

**Download the required database**

The workflow requires a database which can be downloaded and unzipped using these commands.

```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-clone-validation/wf-clone-validation-db.tar.gz
tar -xzvf wf-clone-validation-db.tar.gz
```
The location of the database will need to be specified with the --db_directory parameter. 

## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)