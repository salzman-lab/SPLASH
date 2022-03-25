

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!


Attempt | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Seconds | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269


# stringstats
Motivation

# Prerequisites
1. Install Java.
2. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`).

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


In this example bowtie2 samplesheet, the output anchors and targets will be aligned to the *mm10* and *hg38* indices.
```
/references/mm10/mm10
/references/hg38/hg38
```


## Optional Inputs
*`--anchors_file`*

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
    --anchors_file anchors.txt
    -profile test \
    -r main \
    -latest
```


## Parameters
| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --kmer_size | Length of sequences for anchors and targets | 27 |

`get_anchors`

| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --n_iterations        | Number of chunks of reads to process        | 100 |
| --chunk_size          | Number of reads per chunk to process        | 10000 |
| --target_counts_threshold | Number of unique targets per anchor to initialise phase 1 | 3 |
| --anchor_counts_threshold | Number of total anchor counts to initialise phase 1 | 5 |
| --anchor_freeze_threshold | Maximum number of candidate anchors to store and calculate at a time | 100000 |
| --anchor_score_threshold  | Minimum number of candidate anchors with scores required to calculate anchor significance scores | 100000 |
| --anchor_mode | Mode of fetching candidate anchors from reads, options: `chunk`, `tile` | `tile` |
| --window_slide | If `--anchor_mode tile`, the number of bases to slide across the read to fetch the candidate anchors. If `--window_slide 5`, candidate anchors will start at positions [0,4,9,...] | 5 |
| --use_std | Boolean value if the anchor significance scores should be computed with standard deviation, options: `true`, `false` | `true` |
| --compute_target_distance | Boolean value if the target distance should be computed upon the encounter of a target. If `--compute_target_distance false`, the target distance of a new target will be assigned a conservative estimate of 1, options: `true`, `false` | `false` |

`parse_anchors`

| Argument              | Description       | Default  |
| --------------------- | ----------------- |----------|
| --consensus_length | Maximum length of candidate consensus sequences used to build the final consensus sequence | 200 |
| --direction | The relative direction to search for candidate consensus sequences and targets, options: `up`, `down` | `down` |
| --use_read_length | Boolean value if the looklength should be calculated as a function of read length, options: `true`, `false` | `true` |
| --looklength | If `--use_read_length false`, this is the distance to look for candidate consensus sequences and targets | 0 |
| --num_keep_anchors | Maxiumum number of fastq reads to parse | 4000000 |


## Citations



This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
>

