# Running on your own data
To run this pipeline on your own samplesheet, use a command similar to:
```bash
nextflow run kaitlinchaung/nomad \
    -profile <<singularity/conda/docker>> \
    -r main \
    -latest \
    -resume \
    --input <<YOUR_SAMPLESHEET>> \
    --element_annotations_samplesheet <<ELEMENT_ANN_SAMPLESHEET>>
```
# Inputs
## Required Inputs
*`--input`*

The input samplesheet should be a comma-separated file with no header, consisting of:
1. full paths to gzip fastq files to analyze (required)
2. group IDs of type integer, corresponding to experimental groups of each fastq file (optional)

If paired end sequencing data is being used, please only use files from only Read 1 or files from only Read 2.

In this example samplesheet, 4 fastq files are being analyzed in supervised mode.
```
/data/file1.fastq.gz,-1
/data/file2.fastq.gz,-1
/data/file3.fastq.gz,1
/data/file4.fastq.gz,1
```
In this example samplesheet, 4 fastq files are being analyzed, in unsupervised mode.
```
/data/file1.fastq.gz
/data/file2.fastq.gz
/data/file3.fastq.gz
/data/file4.fastq.gz
```

*`--element_annotations_samplesheet`*

This parameter is a full path to a samplesheet of bowtie2 indices, used in the element annotations step. The default set of bowtie2 indices used in the NOMAD manuscript can be downloaded [here](https://zenodo.org/record/6809531#.YsfR_OzMJTY).

The element annotation samplesheet must not have a header, and it must contain the full path to each bowtie2 index, including the index stem.

Below are general guidlines to creating the element annotation samplesheet:

1. Navigate to a NOMAD index directory.
```
mkdir -p /home/Documents/nomad
cd /home/Documents/nomad
```
2. Download [indices](https://zenodo.org/record/6809531#.YsfR_OzMJTY).
```
wget https://zenodo.org/record/6809531/files/nomad_element_annotation_indices.tar.gz?download=1
```
3. Unpack indices into the index directory.
```
tar -zxvf nomad_element_annotation_indices.tar.gz
```
4. Create the samplesheet, where each line is the full path to each subdirectory from `nomad_element_annotation_indices`, including the reference stems.

For example, if you downloaded `nomad_element_annotation_indices` into `/home/Documents/nomad`,
then your samplesheet would look like the following. Please note that the reference stem--NOT the directory path--is required, otherwise this step will fail.
```
/home/Documents/nomad/nomad_element_annotation_indices/dfam_te_eukaryota/dfam_te_eukaryota
/home/Documents/nomad/nomad_element_annotation_indices/direct_repeats/direct_repeats
/home/Documents/nomad/nomad_element_annotation_indices/escherichia_phage_phiX174/escherichia_phage_phiX174
/home/Documents/nomad/nomad_element_annotation_indices/eukaryota_its1_itstonedb/eukaryota_its1_itstonedb
...
```
5. Pass in the full path to the samplesheet as a parameter of your run.
```
nextflow run kaitlinchaung/stringstats \
    -profile singularity \
    --input /home/data/samplesheet_COVID.csv \
    --element_annotations_samplesheet /home/data/indices_samplesheet.csv \
    -latest
```

Note: Sherlock users who have access to the horence Oak directory do not need to specify this parameter; it will default to a prebuilt-samplesheet on Oak.

## Bypassing significance calculations
*`--anchors_file`*

To input a list of anchors of interest and skip significant calculations, provide this parameter. This can be used if a user wants to generate contingency tables, consensus files, and/or annotation-based outputs from a specific list of anchors.

The anchors file should be a 1 column file with a header, consisting of a list of anchor sequences of interest, with one anchor per line. An example:
```
anchor
AAAAAAAAAA
CCCCCCCCCC
GGGGGGGGGG
```
An example run command with this optional input. Please note that the samplesheet must be provided as well:
```
nextflow run kaitlinchaung/nomad \
    --input samplesheet.csv \
    --anchors_file anchors.txt \
    --element_annotations_samplesheet /home/data/indices_samplesheet.csv \
    -r main \
    -latest
```

# Parameters

Please note that input parameters should be passed with the a double-hypen, while Nextflow-specific parameters should be passed with a single hyphen. Parameters that are not explicitly defined will be set to the defaults below. All below parameters will have a default value assigned to them.

For example:
```
nextflow run kaitlinchaung/nomad
    --input input.txt \
    -r main \
    -latest \
    --run_annotations true
```

## Pipeline Management Parameters
By default, NOMAD performs the following steps:
1. Anchor and target preprocessing (k-mer fetching, counting, abundance filtering)
2. Anchor significance calculations
3. Consensus building
4. Element annotations

The following parameters will add or remove steps from the default pipeline, they should be modified depending on use case. By default, all of these values are set to `false`.

```
--run_trimming <true, false>                Run fastq trimming with TrimGalore.
--run_umitools <true, false>                Run UMITools extraction on 10X fastqs.
--run_decoy <true, false>                   Run abundance control analysis.
--run_annotations <true, false>             Run genome and splicing annotations, and create .heatmaps
--run_anchor_target_counts <true, false>    Create anchor-target counts contingency table.
--run_pvals_only <true, false>              Only run steps 1-2 above.
--run_anchor_heatmaps <true, false>         Create anchor heatmaps, using NOMAD output files.
```

## Anchor Parameters
```
--is_RNAseq <true, false>                   Omit processing of anchors starting with "TTT", for RNAseq-specific analyses
--use_read_length <true, false>             Determine anchor-target distance as a function of read length (default: true)
--lookahead <int>                           Anchor-target distance if `--use_read_length false` (default: null)
--kmer_size                                 Length of anchors and targets (default: 27)
--num_reads_first_pass <int>                Number of FASTQ reads to process on first pass, for anchor significance calculations (default: 4000000)
--num_reads_first_pass <int>                Number of FASTQ reads to process on second pass, for consensus building (default: 4000000)
--anchor_mode <chunk, tile>                 How anchors are fetched, where chunk anchors are adjacent and tile anchors re overlapping (default: tile)
--window_slide <int>                        Size of sliding windows, if `--anchor_mode tile` (default: 5)
--K_num_hashes <int>                        Number of random hashes in significance caluclations (default: 10)
--L_num_random_Cj <int>                     Number of random c_j in significance calculations (default: 50)
--anchor_count_threshold <int>              Minimum number of total counts required to calculate pvalues for an anchor (default: 50)
--anchor_unique_targets_threshold <int>     Minimum number of unique targets to calculate pvalues for an anchor (default: 1)
--anchor_samples_threshold <int>            Minimum number of samples required to calculate pvalues for an anchor (default: 1)
--anchor_sample_counts_threshold <int>      Minimum number of counts across samples required to calculate pvalues for an anchor (default: 5)
--fdr_threshold <float>                     Threshold to call a significant anchor (default: 0.05)
--consensus_length <int>                    Maximum length of candidate consensus sequences used to build the final consensus sequence (default: 200)
```

## Annotation Parameters
These parameters are only used/required if `--run_annotations true`. Otherwise, they do not need to be defined.
```
--genome_index <dir+stem>                   [bowtie2 genome index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) used in genome annotations
--transcriptome_index <dir+stem>            [bowtie2 transcriptome index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) used in genome annotations
--star_index <dir>                          Path to STAR index used in splicing annotations
--gtf <file>                                Path to GTF used in splcing annotations
--gene_bed <file>                           Path to BED file of genes used in genome annotations
```

For `--genome_index` and `--transcriptome_index`, these parameters should be formatted as the directory + reference stem. For example, the parameter would be `--genome_index /home/references/NAME` if `/home/references` contained files such as `NAME.1.bt2`, `NAME.2.bt2`, etc.
