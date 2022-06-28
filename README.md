
## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

# NOMAD
A statistical, reference-free algorithm subsumes myriad problems in genome science and enables novel discovery

## Prerequisites

1. Install Java.
2. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`).
3. Depending on your use case, install [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation), [`docker`](https://www.docker.com/), or [`singularity`](https://sylabs.io/guides/3.5/user-guide/introduction.html). By using the `docker` or `singularity` nextflow profile, the pipeline can be run within the NOMAD docker container (also available on [dockerhub](https://hub.docker.com/repository/docker/kaitlinchaung/stringstats), which contains all the required dependencies.

## Test Run Command
To test this pipeine, use the command below. The `test` profile will launch a pipeline run with a small dataset.
```bash
nextflow run kaitlinchaung/nomad \
    -profile test \
    -r main \
    -latest
```

## Running on your own data
To run this pipeline on your own samplesheet, use a command similar to:
```bash
nextflow run kaitlinchaung/nomad \
    -profile <<singularity/conda/docker>> \
    -r main \
    -latest \
    -resume \
    --input <<YOUR_SAMPLESHEET>>
```

# Inputs
## Required Inputs
*`--input`*

The input samplesheet should be a comma-separated file with no header, consisting of:
1. full paths to gzip fastq files to analyze (required)
2. group IDs of type integer, corresponding to experimental groups of each fastq file (optional)

If paired end sequencing data is being used, please only use files from only Read 1 or files from only Read 2.

In this example samplesheet, 4 fastq files are being analyzed to compare 2 experimental groups.
```
/data/file1.fastq.gz,1
/data/file2.fastq.gz,1
/data/file3.fastq.gz,2
/data/file4.fastq.gz,2
```
In this example samplesheet, 4 fastq files are being analyzed, without comparing any experimental groups.
```
/data/file1.fastq.gz
/data/file2.fastq.gz
/data/file3.fastq.gz
/data/file4.fastq.gz
```

*`--element_annotation_samplesheet`*

The element annotation samplesheet should be 1-column file with no header, consisting of:
1. full paths to [bowtie2 references](bowtie link), including the reference stem

Currently, this paramter defaults to a set of references that is available to Sherlock users.

In this example bowtie2 samplesheet, the output anchors and targets will be aligned to the *mm10* and *hg38* indices.
```
/references/mm10/mm10
/references/hg38/hg38
```


## Optional Inputs
### *`--anchors_file`*

To bypass the `get_anchors` step and input a list of anchors of interest, provide this parameter. Please note that the samplesheet must be provided as well.

The anchors file should be a 1 column file with a header, consisting of a list of anchor sequences of interest, with one anchor per line. An example:
```
anchor
AAAAAAAAAA
CCCCCCCCCC
GGGGGGGGGG
```
An example run command with this optional input:
```
nextflow run kaitlinchaung/nomad \
    --input samplesheet.csv \
    --anchors_file anchors.txt \
    -r main \
    -latest
```


## Parameters

Please note that input parameters should be passed with the a double-hypen, while Nextflow-specific parameters should be passed with a single hyphen. Parameters that are not explicitly defined will be set to the defaults below.

For example:
```
nextflow run kaitlinchaung/nomad \
    --input input.txt \
    -r main \
    -latest \
    --num_lines 2000
```

| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --skip_trimming | Boolean value to indicate if adaptor trimming should be skipped, options: `true`, `false` | `true` |
| --use_read_length | Boolean value to indicate if the distance between anchor and target is a function of read length, options: `true`, `false` | `true` |
| --lookahead | The distance between anchor and target if `--use_read_length true` | 0 |
| --bowtie2_index | Index used for mapping the fastq reads using bowtie2 and extracting the unmapped reads if `--unmapped true` is set | `NA` |
| --run_get_unmapped | Boolean value to indicate if all reads should be used as input, or only the unmapped reads (based on bowtie2 mapping against provided `--bowtie2_index`). If `--run_get_unmapped false`, all reads will be used in this run; if `--unmapped true`, only the unmapped reads will be used in this run; options: `true`, `false`   | `false` |
| --run_decoy | Boolean value to run the decoy version of the pipeline, where the top 1000 most abudnanta anchors are used as pipeline input, options: `true`, `false` | `false` |
| --run_annotations | Boolean value for running genomea nd splicing annotations, options: `true`, `false` | `true` |
| --run_anchor_target_counts | Boolean value for creating a counts table of anchor-targets by sample for visualization, options: `true`, `false` | `false` |
| --run_pvals_only | Boolean value to complete the pipeline after pvalues are completed, options: `true`, `false` | `false` |

*`fetch_anchors`*

| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --num_reads_first_pass | Maximum number of reads to fetch anchors and targets from | 4M |
| --kmer_size | Length of sequences for anchors and targets | 27 |
| --anchor_mode | Mode by which to fetch anchors and target sequences, options: `chunk`, `tile`| `tile` |
| --window_slide | Size of sliding window to fetch anchors, when in `tile` mode | 5 |

*`get_anchors_and_scores`*
| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --anchor_count_threshold | Minimum number of total counts required to calculate a score for an anchor | 50 |
| --K_num_hashes | Number of random hashes | 10 |
| --L_num_random_Cj | Number of random CJ | 50 |
| --fdr_threshold | Pvalue threshold to call a significant anchor | 0.05 |


*`parse_anchors`*

| Argument              | Description       | Default  |
| --------------------- | ----------------- |----------|
| --num_reads_second_pass | Maximum number of reads to build consensus sequences from | 4M |
| --consensus_length | Maximum length of candidate consensus sequences used to build the final consensus sequence | 200 |
| --direction | The relative direction to search for candidate consensus sequences and targets, options: `up`, `down` | `down` |

*Annotation-related parameters*
**If these parameters are not passed in, the pipeline will default to hg38 versions of these files, hosted on Sherlock. For users without Sherlock access, this step may break.**
| Argument              | Description       |
| --------------------- | ----------------- |
| --genome_index | bowtie2 genome index used in `genome_annotations_*` files |
| --transcriptome_index | bowtie2 transcriptome index used in `genome_annotations_*` files |
| --gene_bed | BED file of annotated genes used in `genome_annotations_*` files |
| --exon_starts_bed | BED file of annotated exon start sites used in `genome_annotations_*` files, each entry is from [a, a] where a = exon start |
| --exon_ends_bed | BED file of annotated exon end sites used in `genome_annotations_*` files, each entry is from [a, a] where a = exon end  |
| --star_index | STAR genome index used in splice junction annotations |
| --gtf | GTF file used in splice junction annotations |
| --ann_AS_gtf | alternative splicing annotation file |


## Citations


This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
>

