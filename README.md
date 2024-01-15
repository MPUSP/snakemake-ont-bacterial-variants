# Snakemake workflow: snakemake-ont-bacterial-variants

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for the identification of variants in bacterial genomes using nanopore long-read sequencing.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Workflow overview

## Installation

### Snakemake

Step 1: Install snakemake with `conda` in a new conda environment.

```
conda create -n <ENV> snakemake
```

Step 2: Activate conda environment with snakemake

```
source /path/to/conda/bin/activate
conda activate <ENV>
```

### Additional tools

**Important note:**

All other dependencies for the workflow are **automatically pulled as `conda` environments** by snakemake, when running the workflow with the `--use-conda` parameter (recommended).
When run without automatically built `conda` environments, all packages need to be installed manually.

## Running the workflow

## Output

## Authors

## References
