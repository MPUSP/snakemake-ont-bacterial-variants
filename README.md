# Snakemake workflow: snakemake-ont-bacterial-variants

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/MPUSP/snakemake-ont-bacterial-variants/workflows/Tests/badge.svg?branch=main)](https://github.com/MPUSP/snakemake-ont-bacterial-variants/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow for the identification of variants in bacterial genomes using nanopore long-read sequencing.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=MPUSP/snakemake-ont-bacterial-variants).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Workflow overview

This workflow provides a simple and easy-to-use framework for the identification of structural and small nucleotide variants in bacterial genomes using nanopore long-read sequencing data. 
The `snakemake-ont-bacterial-variants` workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Quality check of sequencing data (`FastQC`)
2. Filtering of input sequencing data by read length and quality (`Filtlong`)
3. Mapping to reference genome (`NGMLR`)
4. Calling of structural and single nucleotide variants (SVs: `cuteSV` and `Sniffles2`; SNVs: `Clair3` and `Medaka`)
5. Filtering of identified variants (e.g., by variant quality or genomic regions; `BCFtools` and `VCFtools`)
6. Generate report with final results (`R markdown`, `igv-reports`, and `MultiQC`) 

If you would like to contribute, report issues, or suggest features, please get in touch on [GitHub](https://github.com/MPUSP/snakemake-ont-bacterial-variants).

## Installation

### Snakemake

Step 1: Install snakemake with `conda` in a new conda environment.

```
conda create -n <ENV> snakemake
```

Step 2: Activate conda environment with snakemake

```
conda activate <ENV>
```

### Additional tools

**Important note:**

All other dependencies for the workflow are **automatically pulled as `conda` environments** by snakemake, when running the workflow with the `--use-conda` parameter (recommended).
When run without automatically built `conda` environments, all packages need to be installed manually:
- `NanoPlot`
- `MultiQC`
- `Filtlong`
- `NGMLR`
- `minimap2`
- `samtools`
- `bedtools`
- `Medaka`
- `Clair3`
- `cuteSV`
- `Sniffles2`
- `bcftools`
- `VCFtools`
- `r-tidyverse`
- `r-rmarkdown`
- `r-dt`
- `igv-reports`

## Running the workflow

### Input data

The workflow requires the following files to be located in the `data` directory:
1. Whole-genome sequencing data in `*.fastq.gz` format in `data/fastq`
2. Reference genome(s) in `*.fa` format in `data/reference`

Optionally, users can provide:
- Reference genome annotation in `*.gff` format in `data/annotation` (for feature annotation in IGV report)
- A `*.bed` file with genomic regions to ignore for variant calling in `data/masked_region`

Please ensure that the chromosome names in `*.gff` and `*.bed` files are the same as in the reference genome.

Input data files are provided in the `samples.tsv` table, whose location is inidcated in the `config.yml` file. The samplesheet must comply with the following structure:
- `sample` defines the sample name that will be used throughout the workflow and thus needs to be **unique**.
- `fastq` provides the path to the sample's `*.fastq.gz` file.
- `reference` provides the path to the reference genome `*.fa` file (may be the same for several / all samples).
- `annotation` provides the path to the optional reference genome annotation in `*.gff` file (may be the same for several / all samples). If no annotation is provided, you **must** enter `n/a`!
- `masked_regions` provides the path to the optional `*.bed` file for filtering genomic regions (may be the same for several / all samples). If no `*.bed` file is provided, you **must** enter `n/a`!

| sample    | fastq                        | reference                | annotation                  | masked_regions                   |
| --------- | ---------------------------- | ------------------------ | --------------------------- | -------------------------------- |
| <sample1> | data/fastq/<fastq1>.fastq.gz | data/reference/<ref1>.fa | data/annotation/<anno1>.gff | data/masked_region/<region1>.bed |
| <sample2> | data/fastq/<fastq2>.fastq.gz | data/reference/<ref2>.fa | data/annotation/<anno2>.gff | data/masked_region/<region2>.bed |
| ...       | ...                          | ...                      | ...                         | ...                              |
| <sampleN> | data/fastq/<fastqN>.fastq.gz | data/reference/<refN>.fa | data/annotation/<annoN>.gff | data/masked_region/<regionN>.bed |

### Configuration and parameters

Before executing the workflow, you may want to adjust several options and parameters in the default config file `config/config.yml`:
1. Directories:
 * `indir`: Input directory for all input files, `data` by default (see above)
 * `outdir`: Output directory (relative to working directory), `results` by default
2. Sample information:
 * `samples`: Path to samplesheet (relative to working directory), `samplesheet/samples.tsv` by default
 * `libprepkit`: Kit from ONT used for library preparation, e.g. `SQK-NBD114.24`
 * `basecalling_model`: Model used for basecalling of raw sequencing data (required for variant calling using `Medaka`), currently supported models are:
  * `r1041_e82_400bps_sup_v4.2.0`
  * `r1041_e82_400bps_sup_v4.3.0`
3. Tool parameters:
 * The number of cores can be adjusted here for the following tools: `NGMLR`, `NanoPlot`, `MultiQC`, `Medaka`, `Clair3`, `Sniffles2`, and `cuteSV`
 * You may further adjust the run parameters for the following tools (please refer to the reference provided for more details on run parameters):
  * `Filtlong`: By default, reads are filtered for a minimum length of 500 bp and a mean accuracy of at least 90% (Q10), with 90% of the longest and highest-quailty reads to be kept.
  * `Clair3`: Variants are called on all contigs in a haploid-sensitive, ONT-specific mode using `--include_all_ctgs --haploid_sensitive --platform ont`.
  * `cuteSV`: Variants are called with the suggested parameters for ONT data (`--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3`) and the genotyping option enabled (`--genotype`). 
4. Filtering of variants:
 * The variant quality thresholds can be adjusted here for all four variant callers
 * `remove_common_variants`: If `True`, variants which have been identified in all samples with the same reference genome by one tool are filtered out. This is helpful in case all samples derive from a strain, whose genome sequence already differs from the used reference sequence. If `False`, all variants are reported.
5. Reporting options:
 * `igv_region_length`: Neighboring variants with a maximum bp distance indicated here [1 by default] are reported in one region in the IGV variant report. Increasing this parameter will reduce the file size of the resulting IGV HTML report, if hotspots / regions with many variants exist in a sample.

### Execution

To run the workflow from the command line, change the working directory:

```
cd /path/to/snakemake-ont-bacterial-variants
```

Before running the entire workflow, you may perform a dry run using:

```
snakemake --dry-run
```

To run the complete workflow with your files using **`conda`**, execute the following command (the definition of the number of cores is mandatory):

```
snakemake --cores 10 --use-conda
```

You may also run the workflow on the provided test data using:

```
snakemake --cores 10 --use-conda --all-temp --directory .test
```

## Output

### Main output

The most important output files and reports are found in `variant_reports` for each sample in a separate sub-directory `variant_reports/<sample>/`:
- `<sample>_overview.html`: Custom overview HTML report comprising read and mapping statistics as well as identified variants.
- `<sample>_IGV.html`: HTML report with alignment data for all variants listed in `<sample>_overview.html` (generated with `igv-reports`).
- `<sample>.<tool>.vcf`: File with all variants for each tool used (after filtering).

Further, `variants_reports` contains `MultiQC.html`, an interactive HTML report generated by `MultiQC` with read statistics on raw, filtered and aligned reads.

Log files from the generation of above reports can be found in `variant_reports/logs/`:
- `<sample>.collect_vcfs.log`: Copying of variant files to above destinations
- `<sample>.igv_reports.log`: Generation of IGV HTML report
- `<sample>.prepare_vcfs.log`: Merging of neighboring variants in genomic regions for `igv-reports`
- `<sample>.report.log`: Generation of final overview HTML report using R Markdown

### Additional output files

In addition, the workflow generates the following output files:

<details markdown="1">
<summary>filtered_reads</summary>

- `<sample>.fastq.gz`: Length- and quality-filtered reads
- `logs/<sample>.log`: Filtering log file(s)

</details>

<details markdown="1">
<summary>mapping</summary>

- `<sample>.bam`: Alignment file of filtered reads to sample-specific reference genome (produced by `NGMLR`)
- `<sample>.bam.bai`: Alignment index file
- `logs/<sample>.*.log`: Log files for `NGMLR`, `samtools`, and `bedtools` (calculation of genome coverage and read end distributions [as temporary files only])

</details>

<details markdown="1">
<summary>qc</summary>

Read statistics on all raw, filtered and aligned reads are found here (generated with `NanoPlot`):
- `raw_reads/<sample>/`: Directory with output files with read statistics on all **raw reads** (log files in `filtered_reads/logs/<sample>.log`).
- `filtered_reads/<sample>/`: Directory with output files with read statistics on all **filtered reads** (log files in `filtered_reads/logs/<sample>.log`).
- `aligned_reads/<sample>/`: Directory with output files with read statistics on all **aligned reads** (log files in `aligned_reads/logs/<sample>.log`).

The raw output from `MultiQC` - based on above `NanoPlot` results - is located in the directory `multiqc`, containing the interactive HTML report `multiqc/multiqc_report.html`.

</details>

<details markdown="1">
<summary>SNV</summary>

For both tools (`medaka` and `clair3`):
- `<tool>/<sample>/`: Directory with all **raw** output files from variant calling tool
- `<tool>/<sample>.vcf`: File with identified variants
- `<tool>/<sample>.filtered.vcf`: File with variants filtered by variant quality, overlap with other samples from same reference genome (only if `remove_common_variants: True` in `config.yml`), and by genomic region if provided (compare `masked_region` in samplesheet)
- `<tool>/logs/<sample>.*`: Standard error (`stderr`) and standard output (`stdout`) files for the identification of variants using either `medaka` or `clair3` 
- `<tool>/logs/<sample>_filtering.log`: Log file produced by `VCFtools` for final filtering step to generate `<tool>/<sample>.filtered.vcf`

If `remove_common_variants: True` in `config.yml`, the following files are additionally produced:
- `<tool>/common_variants.vcf`: File with variants called by the tool which are shared by all samples with the same reference genome
- `<tool>/common_variants.txt`: Tab-delimited file deduced from `<tool>/common_variants.vcf` listing contig and variant start position only
- `<tool>/logs/<sample>_indexing.*`: Standard error (`stderr`) and standard output (`stdout`) for indexing of `*.vcf` files using `bcftools`
- `<tool>/logs/bcftools.*`: Standard error (`stderr`) and standard output (`stdout`) for identification of shared variants across samples uisng `bcftools`

For `clair3`, data from the downloaded model is additionally found in the directory `clair3/model/`.

</details>

<details markdown="1">
<summary>SV</summary>

For both tools (`cutesv` and `sniffles2`):
- `<tool>/<sample>.vcf`: File with identified variants
- `<tool>/<sample>.filtered.vcf`: File with variants filtered by variant quality, overlap with other samples from same reference genome (only if `remove_common_variants: True` in `config.yml`), and by genomic region if provided (compare `masked_region` in samplesheet)
- `<tool>/logs/<sample>.log`: Log file for the identification of variants using either `cutesv` or `sniffles2`
- `<tool>/logs/<sample>_filtering.log`: Log file produced by `VCFtools` for final filtering step to generate `<tool>/<sample>.filtered.vcf`

If `remove_common_variants: True` in `config.yml`, the following files are additionally produced:
- `<tool>/common_variants.vcf`: File with variants called by the tool which are shared by all samples with the same reference genome
- `<tool>/common_variants.txt`: Tab-delimited file deduced from `<tool>/common_variants.vcf` listing contig and variant start position only
- `<tool>/logs/<sample>_indexing.*`: Standard error (`stderr`) and standard output (`stdout`) for indexing of `*.vcf` files using `bcftools`
- `<tool>/logs/bcftools.*`: Standard error (`stderr`) and standard output (`stdout`) for identification of shared variants across samples uisng `bcftools`

For `cutesv`, all **raw** output files from the initial variant calling can be found in the directory `cutesv/<sample>/`.

</details>

## Authors

Dr. Thomas Fabian Wulff
- Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
- ORCID: https://orcid.org/0000-0001-7166-0899

Please also visit the MPUSP GitHub page at https://github.com/MPUSP for further information on this workflow and other projects!

## References

The following software and tools have been used in this workflow:

<details markdown="1">
<summary>BCFtools</summary>

> Danecek, P., Bonfield, J.K., Liddle, J. et al. _Twelve years of SAMtools and BCFtools_. GigaScience 10(2), giab008, 2021. (https://doi.org/10.1093/gigascience/giab008)

</details>

<details markdown="1">
<summary>bedtools</summary>

> Quinlan, A.R., Hall, I.M., _BEDTools: a flexible suite of utilities for comparing genomic features_. Bioinformatics 26(6), 841-842, 2010. ( https://doi.org/10.1093/bioinformatics/btq033)

</details>

<details markdown="1">
<summary>Clair3</summary>

> Zheng, Z., Li, S., Su, J. et al. _Symphonizing pileup and full-alignment for deep learning-based long-read variant calling_. Nature Computational Science 2, 797–803, 2022. (https://doi.org/10.1038/s43588-022-00387-x)

</details>

<details markdown="1">
<summary>cuteSV</summary>

> Jiang, T., Liu, Y., Jiang, Y. et al. _Long-read-based human genomic structural variation detection with cuteSV_. Genome Biology 21, 189, 2020. (https://doi.org/10.1186/s13059-020-02107-y)

</details>

<details markdown="1">
<summary>Filtlong</summary>

> https://github.com/rrwick/Filtlong

</details>

<details markdown="1">
<summary>igv-reports</summary>

> https://github.com/igvteam/igv-reports

</details>

<details markdown="1">
<summary>Medaka</summary>

> https://github.com/nanoporetech/medaka

</details>

<details markdown="1">
<summary>MultiQC</summary>

> Ewels, P., Magnusson, M., Lundin, S., et al. _MultiQC: summarize analysis results for multiple tools and samples in a single report_. Bioinformatics 32(19) 3047–3048, 2016. (https://doi.org/10.1093/bioinformatics/btw354)

</details>

<details markdown="1">
<summary>NanoPlot</summary>

> Coster, W.D., Rademakers, R. _NanoPack2: population-scale evaluation of long-read sequencing data_. Bioinformatics 39(5), btad311, 2023. (https://doi.org/10.1093/bioinformatics/btad311)

</details>

<details markdown="1">
<summary>NGMLR</summary>

> Sedlazeck, F.J., Rescheneder, P., Smolka, M. et al. _Accurate detection of complex structural variations using single-molecule sequencing_. Nature Methods 15, 461–468, 2018. (https://doi.org/10.1038/s41592-018-0001-7)

</details>

<details markdown="1">
<summary>SAMtools</summary>

> Li, H., Handsaker, B., Wysoker, A. et al. _The Sequence Alignment/Map format and SAMtools_. Bioinformatics 25(16), 2078–2079, 2009. (https://doi.org/10.1093/bioinformatics/btp352)

</details>

<details markdown="1">
<summary>Snakemake</summary>

> Mölder, F., Jablonski, K.P., Letcher, B. et al. _Sustainable data analysis with Snakemake_. F1000Research 10:33, 2021. (https://doi.org/10.12688/f1000research.29032.2)

</details>

<details markdown="1">
<summary>Sniffles2</summary>

> Smolka, M., Paulin, L.F., Grochowski, C.M. et al. _Detection of mosaic and population-level structural variants with Sniffles2_. Nature Biotechnology 2024. (https://doi.org/10.1038/s41587-023-02024-y)

</details>

<details markdown="1">
<summary>VCFtools</summary>

> Danecek, P., Auton, A., Abecasis, G. et al. _The variant call format and VCFtools_. Bioinformatics 27(15), 2156–2158, 2011. (https://doi.org/10.1093/bioinformatics/btr330)

</details>