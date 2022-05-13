test

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

test
# stringstats
Motivation

# Prerequisites
**For sherlock users, more details found at https://github.com/kaitlinchaung/stringstats/blob/main/nextflow.md**

1. Install Java.
2. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`).
3. Install [`docker`](https://www.docker.com/) or [`singularity`](https://sylabs.io/guides/3.5/user-guide/introduction.html). By using the `docker` or `singularity` nextflow profile, the pipeline can be run within the stringstats docker container (also available on [dockerhub](https://hub.docker.com/repository/docker/mariekevromman/stringstats)), which contains all the required dependencies.

# Test Run Command
To test this pipeine, use the command below. The `test` profie will launch a pipeline run with a small dataset.
```bash
nextflow run kaitlinchaung/stringstats \
    -profile test \
    -r main \
    -latest
```

# Inputs
## Required Inputs
*`--input`*

The input samplesheet should be a comma-separated file with no header, consisting of:
1. full paths to gzip fastq files to analyze (required)
2. group IDs of type integer, corresponding to experimental groups of each fastq file (optional)

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

*`--bowtie2_samplesheet`*

The bowtie2 samplesheet should be 1-column file with no header, consisting of:
1. full paths to [bowtie2 references](bowtie link), including the reference stem

Currently, this paramter defaults to a set of references that is available to Sherlock users.

In this example bowtie2 samplesheet, the output anchors and targets will be aligned to the *mm10* and *hg38* indices.
```
/references/mm10/mm10
/references/hg38/hg38
```


## Optional Inputs
### *`--anchors_file`*

To bypass the `get_anchors` step and input a list of anchors of interest, provide this parameter.

The anchors file should be a 1 column file with a header, consisting of a list of anchor sequences of interest, with one anchor per line. An example:
```
anchor
AAAAAAAAAA
CCCCCCCCCC
GGGGGGGGGG
```
An example run command with this optional input:
```
nextflow run kaitlinchaung/stringstats \
    --input samplesheet.csv \
    --anchors_file anchors.txt \
    -r main \
    -latest
```


## Parameters

Please note that input parameters should be passed with the a double-hypen, while Nextflow-specific parameters should be passed with a single hyphen. Parameters that are not explicitly defined will be set to the defaults below.

For example:
```
nextflow run kaitlinchaung/stringstats \
    --input input.txt \
    -r main \
    -latest \
    --num_lines 2000
```

| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --use_read_length | Boolean value to indicate if the distance between anchor and target is a function of read length, options: `true`, `false` | `true` |
| --lookahead | The distance between anchor and target if `--use_read_length true` | 0 |
| --unmapped | Boolean value to indicate if all reads should be used as input, or only the unmapped reads (based on bowtie2 mapping against provided `--bowtie2_index`). If `--unmapped false`, all reads will be used in this run; if `--unmapped true`, only the unmapped reads will be used in this run; options: `true`, `false`   | `false` |
| --bowtie2_index | Index used for mapping the fastq reads using bowtie2 and extracting the unmapped reads if `--unmapped true` is set | `NA` |

*`fetch_anchors`*

| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --num_lines | Maximum number of reads to fetch anchors and targets from | no maximum aka 0 |
| --kmer_size | Length of sequences for anchors and targets | 27 |
| --anchor_mode | Mode by which to fetch anchors and target sequences, options: `chunk`, `tile`| `tile` |
| --window_slide | Size of sliding window to fetch anchors, when in `tile` mode | 5 |

*`get_anchors_and_scores`*
| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --anchor_count_threshold | Minimum number of total counts required to calculate a score for an anchor | 50 |
| --distance_type | options: `hamming`, `lev` | `lev` |
| --max_targets | Maximum number of targets per anchor for score calculation | 50 |
| --max_dist | Maximum distance allowed between targets in score calculation | 10 |
| --bonfer | Number of corrections | 10 |
| --pval_threshold | Pvalue threshold to call a significant anchor | 0.1 |


*`parse_anchors`*

| Argument              | Description       | Default  |
| --------------------- | ----------------- |----------|
| --num_parse_anchors_reads | Maximum length of candidate consensus sequences used to build the final consensus sequence | 4000000 |
| --consensus_length | Maximum length of candidate consensus sequences used to build the final consensus sequence | 200 |
| --direction | The relative direction to search for candidate consensus sequences and targets, options: `up`, `down` | `down` |

*Annotation-related parameters*
| Argument              | Description       |
| --------------------- | ----------------- |
| --genome_index | bowtie2 genome index used in `genome_annotations_*` files |
| --transcriptome_index | bowtie2 transcriptome index used in `genome_annotations_*` files |
| --gene_bed | BED file of annotated genes used in `genome_annotations_*` files |
| --exon_starts_bed | BED file of annotated exon start sites used in `genome_annotations_*` files |
| --exon_ends_bed | BED file of annotated exon end sites used in `genome_annotations_*` files |
| --star_index | STAR genome index used in splice junction annotations |
| --star_index | GTF file used in splice junction annotations |


## Citations


This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
>

