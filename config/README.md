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

| sample     | fastq                         | reference                 | annotation                   | masked_regions                    |
| :--------- | :---------------------------- | :------------------------ | :--------------------------- | :-------------------------------- |
| \<sample1> | data/fastq/\<fastq1>.fastq.gz | data/reference/\<ref1>.fa | data/annotation/\<anno1>.gff | data/masked_region/\<region1>.bed |
| \<sample2> | data/fastq/\<fastq2>.fastq.gz | data/reference/\<ref2>.fa | data/annotation/\<anno2>.gff | data/masked_region/\<region2>.bed |
| ...        | ...                           | ...                       | ...                          | ...                               |
| \<sampleN> | data/fastq/\<fastqN>.fastq.gz | data/reference/\<refN>.fa | data/annotation/\<annoN>.gff | data/masked_region/\<regionN>.bed |

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