

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!


# stringstats
Motivation

# Prerequisites
1. Install Java.
2. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`).
3. Install [`docker`](https://www.docker.com/) or [`singularity`](https://sylabs.io/guides/3.5/user-guide/introduction.html). Sherlock has alread preinstalled Singularity.
# Test Run Command
To test this pipeine, use the command below. The `test_consensus` profile will launch a pipeline run with a small dataset hosted on Sherlock.
```bash
nextflow run kaitlinchaung/stringstats \
    -profile test_consensus,singularity,sherlock \
    -r consensus \
    -latest \
    -resume
```
# Run your own data
To run your own data with a looklength=24, you would use something like this below command.
To pass in parameters that are different than the default, use double hyphens.
```bash
nextflow run kaitlinchaung/stringstats \
    -profile singularity,sherlock \
    -r consensus \
    -latest \
    -resume \
    --input samplesheet.csv \
    --anchors.tsv anchors.tsv \
    --looklength 24
```
# Inputs
## Required Inputs
*`--input`*

The input samplesheet should be a comma-separated file with no header, consisting of:
    1. full paths to gzip fastq files to analyze, path basenames must be unique ie file1, file2, file3

This file must have an extension of `.csv`.

In this example samplesheet, 4 fastq files are being analyzed.
```
/data/file1.fastq.gz
/data/file2.fastq.gz
/data/file3.fastq.gz
/data/file4.fastq.gz
```

*`--anchors_file`

To bypass the `get_anchors` step and input a list of anchors of interest, provide this parameter.

The anchors file should be a 1 column file with a header, consisting of a list of anchor sequences of interest, with one anchor per line. An example:
```
anchor
AAAAAAAAAA
CCCCCCCCCC
GGGGGGGGGG
```


## Parameters

Please note that input parameters should be passed with the a double-hypen, while Nextflow-specific parameters should be passed with a single hyphen. Parameters that are not explicitly defined will be set to the defaults below.

For example:
```
nextflow run kaitlinchaung/stringstats \
    --input input.txt \
    -r main \
    -latest \
    --n_itrations 200 \
    --chunk_size 50000
```

| Argument              | Description       | Default  |
| --------------------- | ----------------- |--------- |
| --kmer_size | Length of sequences for anchors and targets | 27 |
| --consensus_length | Maximum length of candidate consensus sequences used to build the final consensus sequence | 200 |
| --direction | The relative direction to search for candidate consensus sequences and targets, options: `up`, `down` | `down` |
| --use_read_length | Boolean value if the looklength should be calculated as a function of read length, options: `true`, `false` | `true` |
| --looklength | If `--use_read_length false`, this is the distance to look for candidate consensus sequences and targets | 0 |
| --num_keep_anchors | Maxiumum number of fastq reads to parse | 4000000 |


# Outputs


`parse_anchors`
1. `consensus_anchors/`
    1. `*.fasta`
        * Fasta file of the consensus sequences of each significant anchor sequence
    2. `*_counts.tab`
        * For each consensus sequence, the frequency count of the consensus base at each position
    3. `*_fractions.tab`
        * For each consensus sequence, the fraction of occurrence of the consensus base at each position
2. `target_counts/`
    * For each fastq file, a counts file for each significant anchor and all of their valid targets




## Citations



This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
>

