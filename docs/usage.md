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

The input samplesheet should be a comma-separated file with no header.

If paired end sequencing data is being used, please only use files from only Read 1 or files from only Read 2.

Here is an example samplesheet for bulk RNAseq or SS2 fastq files. The first column is the sample--this is where metadata information can be incorporated. The second column is the unique fastq file identifier. The third column is the full path to the fastq.
```
muscle_cells,sample_1,/data/sample_1.fastq.gz
muscle_cells,sample_2,/data/sample_2.fastq.gz
brain_cells,sample_3,/data/sample_3.fastq.gz
brain_cells,sample_4,/data/sample_4.fastq.gz
```

Here is an example samplesheet for 10X fastq files. The first column is sample or channel name, where channels can be used to process many fastq files as part of the same group. The second column is the R1 fastq file, and the third column is the paired R2 fastq file. In this case, `--cells_split_across_lanes true` because many FASTQ files belong to the same channel.
```
muscle,/data/muscle_L001_R1.fastq.gz,/data/muscle_L001_R2.fastq.gz
muscle,/data/muscle_L002_R1.fastq.gz,/data/muscle_L002_R2.fastq.gz
muscle,/data/muscle_L003_R1.fastq.gz,/data/muscle_L003_R2.fastq.gz
muscle,/data/muscle_L004_R1.fastq.gz,/data/muscle_L004_R2.fastq.gz
brain,/data/brain_L001_R1.fastq.gz,/data/brain_L001_R2.fastq.gz
brain,/data/brain_L002_R1.fastq.gz,/data/brain_L002_R2.fastq.gz
brain,/data/brain_L003_R1.fastq.gz,/data/brain_L003_R2.fastq.gz
brain,/data/brain_L004_R1.fastq.gz,/data/brain_L004_R2.fastq.gz
```

Here is an example samplesheet for 10X fastq files. The first column is sample or channel name, where channels can be used to process many fastq files as part of the same group. The second column is the R1 fastq file, and the third column is the paired R2 fastq file. In this case, `--cells_split_across_lanes false` because one FASTQ file correspondds to one channel.
```
muscle,/data/muscle_R1.fastq.gz,/data/muscle_R2.fastq.gz
brain,/data/brain_R1.fastq.gz,/data/brain_R2.fastq.gz
liver,/data/liver_R1.fastq.gz,/data/liver_R2.fastq.gz
lung,/data/lung_R1.fastq.gz,/data/lung_R2.fastq.gz
```

10X runs also require another parameter, `--cell_barcode_samplesheet`. The first column must contain the metadata information, such as cell type. The second column must contain the cell barcode, and the third column must contain the channel. Here is an example of a `--cell_barcode_samplesheet`.
```
capillary,AAAA,muscle
capillary,CCCC,muscle
macrophage,GGGG,muscle
macrophage,TTTT,muscle
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

# Parameters

Please note that input parameters should be passed with the a double-hypen, while Nextflow-specific parameters should be passed with a single hyphen. Parameters that are not explicitly defined will be set to the defaults below. All below parameters will have a default value assigned to them.

For example:
```
nextflow run kaitlinchaung/nomad
    --input input.txt \
    -r main \
    -latest \
    --is_10X true
```

## Pipeline Management Parameters
By default, NOMAD performs the following steps:
1. Anchor and target preprocessing (k-mer fetching, counting, abundance filtering)
2. Anchor significance calculations
3. Element annotations
4. Genome annotations
5. Summary of significant anchors and their annotations

The following parameters will add or remove steps from the default pipeline, they should be modified depending on use case. By default, all of these values are set to `false`.

```
--run_umitools <true, false>                Run UMITools extraction on 10X fastqs.
--run_control <true, false>                   Run abundance control analysis.
--run_pvals_only <true, false>              Only run steps 1-2 above.
```

## Parallelization Parameters
```
--stratify_level <int>                      The number of parallelizations to utilize during p value computations, will be equal to 4^stratify_level parallelizations (default: 3)
--anchor_batch_size <int>                   The number of anchors used in p value computation per batch, lower values will use less memory (default: 1000)
```

## Anchor Parameters
```
--run_unsupervised_pvals <true, false>      If pvalues should be calculated without metadata information (default: false)
--is_RNAseq <true, false>                   Omit processing of anchors starting with "TTT", for RNAseq-specific analyses (default: false)
--use_read_length <true, false>             Determine anchor-target distance as a function of read length (default: true)
--lookahead <int>                           Anchor-target distance if `--use_read_length false` (default: null)
--kmer_size                                 Length of anchors and targets (default: 27)
--num_reads_first_pass <int>                Number of FASTQ reads to process on first pass, for anchor significance calculations (default: 4000000)
--anchor_mode <chunk, tile>                 How anchors are fetched, where chunk anchors are adjacent and tile anchors re overlapping (default: tile)
--window_slide <int>                        Size of sliding windows, if `--anchor_mode tile` (default: 5)
--K_num_hashes <int>                        Number of random hashes in significance caluclations (default: 10)
--L_num_random_Cj <int>                     Number of random c_j in significance calculations (default: 50)
--anchor_count_threshold <int>              Minimum number of total counts required to calculate pvalues for an anchor (default: 50)
--anchor_unique_targets_threshold <int>     Minimum number of unique targets to calculate pvalues for an anchor (default: 1)
--anchor_samples_threshold <int>            Minimum number of samples required to calculate pvalues for an anchor (default: 1)
--anchor_sample_counts_threshold <int>      Minimum number of counts across samples required to calculate pvalues for an anchor (default: 5)
--fdr_threshold <float>                     Threshold to call a significant anchor (default: 0.05)
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

## 10X Parameters
```
--is_10X <true, false>                      Must be `true` to run workflow in 10X mode (default: false)
--run_umitools <true, false>                Run umitools whitelist and extract on fastq files (default: true)
--cells_split_across_lanes <true, false>    If the same cells have been sequenced across many sequencing lanes
```

For `--genome_index` and `--transcriptome_index`, these parameters should be formatted as the directory + reference stem. For example, the parameter would be `--genome_index /home/references/NAME` if `/home/references` contained files such as `NAME.1.bt2`, `NAME.2.bt2`, etc.
