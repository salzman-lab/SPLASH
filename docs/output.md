# kaitlinchaung/nomad: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.


## File Output Descriptors

### `./anchors_pvals.tsv`
1. anchor: the anchor in question
2. pv_hand: uncorrected p values with handcrafted di's
3. pv_hash: uncorrected p values from random hash functions
4. pv_hand_sheetCjs: uncorrected p values from handcrafted di’s using the samplesheet 5. c_j’s
5. pv_hash_sheetCjs: uncorrected p values from random hash functions using the samplesheet c_j’s
6. effectSize_randCjs: measure of effect size listed in previous section, using random cjs and corresponding maximizing random hash (paired with pv_hash)
7. effectSize_sheetCjs: measure of effect size above, using sheet cjs and corresponding maximizing random hash (analogous to pv_hash_sheetCjs)
8. OptHash: index of optimizing hash for pv_hash computation, mmh3 hash library with this seed
9. M: total number of reads used in the contingency table
10. entropy: entropy of total target distribution for this anchor (a label agnostic measure)
11. mu_ham: mu for handcrafted di's: for each target sequence, we take di to be the target's hamming distance to the most abundant target. This yields the handcrafted d_i’s, which take values between 0 and 27 (the k-mer size). For mu_ham, we take the average of all of these over this given anchor. This measure indicates diversity of target sequences.
12. mu_lev: same as mu_ham, but with levenshtein distance instead of hamming.
13. pv_hash_both: 2*min(pv_hash, pv_hash_sheetCjs), Bonferonni corrected p value
14. pv_hash_both_corrected: BY corrected q values for pv_hash_both across all anchors for this dataset (Benjamini Yekutieli, BH with correction for arbitrary dependence)
cj_rand_${sample}: optimizing cj for pv_hash for the given sample

### `./element_annotations/annotated_[anchors/targets].tsv`
1. [anchor/target]: the sequence in question
2. For all input references:
    1. [anchor/target]\_hits\_${reference}: the alignment of the sequence to the reference
    2. [anchor/target]\_hits\_pos_${reference}: the alignment position of the sequence to the reference

### `./genome_annotations/genome_annotations_[anchor/target/consensus_fastq_id].tsv`
1. [anchor/target/anchor of the consensus]: the sequence in question
2. local_strand: when aligned in bowtie2 local mode, the strand of the sequence genome alignment
3. local_gene: when aligned in bowtie2 local mode, a unique list of annotated genes that have any intersection with the sequence genome alignment
4. local_gene_MAPQ:  when aligned in bowtie2 local mode, the MAPQ score of the sequence genome alignment
5. local_transcript: when aligned in bowtie2 local mode, the sequence alignment to the reference transcriptome
6. local_transcriptome_MAPQ: when aligned in bowtie2 local mode, the MAPQ score of the sequence transcriptome alignment
7. end_to_end_strand: when aligned in bowtie2 end-to-end mode, the strand of the sequence genome alignment
8. end_to_end_gene: when aligned in bowtie2 end-to-end mode, a unique list of annotated genes that have any intersection with the sequence genome alignment
9. end_to_end_gene_MAPQ: when aligned in bowtie2 end-to-end mode, the MAPQ score of the sequence genome alignment
10. end_to_end_transcript: when aligned in bowtie2 end-to-end mode, the sequence alignment to the reference transcriptome
11. end_to_end_transcriptome_MAPQ: when aligned in bowtie2 end-to-end mode, the MAPQ score of the sequence transcriptome alignment

### `./consensus_anchors/splicing_annotations/consensus_called_exons.tsv`
1. called_exon_chr: the chromosome of the called exon boundary
2. called_exon_id: if the called exon boundary represents the start or end position of the called exon
3. called_exon_id_position: the genomic position of the called exon boundary
4. ann_exon_gene: if the called exon boundary is within report 1base pair of an annotated exon boundary, report tha annotated gene name
5. ann_exon_strand: if the called exon boundary is within report 1base pair of an annotated exon boundary, report that annotated gene strand
6. sample: sample ID, corresponding to the fastq file from which the consensus sequence was built
7. anchor: anchor sequence of the consensus sequence
8. consensus: consensus sequence from which the called exons are called from
9. anchor_local_gene: when the anchor is aligned in bowtie2 end-to-end mode, a unique list of annotated genes that have any intersection with the anchor genome alignment
10. anchor_end_to_end_gene: when the anchor is aligned in bowtie2 local mode, a unique list of annotated genes that have any intersection with the anchor genome alignment
11. consensus_gene: a list of unique genes that have any intersection with the consensus` called exons
12. consensus_reported_aligment: the consensus sequence’s alignment, as reported by STAR

### `./summary.tsv`
1. anchor: anchor sequence
2. target: target sequence
3. pv_hash
4. pv_hand
5. effect_size_randCjs
6. optHash
7. M
8. entropy mu_ham
9. mu_lev
10. pv_hash_both
11. pv_hash_both_corrected
12. rcAnchor: reverse complement of the anchor sequence
13. rcTarget: reverse complement of the target sequence
14. anchor_matches_rcTarget: if the anchor has any match to any reverse complemented target
15. rcAnchor_matches_target: if the target has any match to any reverse complemented anchor
16. anchor_top_ann: the highest priority element annotation for the anchor
17. anchor_top_ann_hit: the alignment of the highest priority element annotation for the anchor
18. anchor_top_ann_hit_pos: the position of the highest priority element annotation for the anchor
19. anchor_num_ann: the number of element annotations that the anchor has
20. anchor_annotation_source: if the above anchor alignments were performed with bowtie2 or blast
21. target_top_ann  target_top_ann_hit: the highest priority element annotation for the target
22. target_top_ann_hit_pos: the alignment of the highest priority element annotation for the target
23. target_num_ann: the number of element annotations that the target has
24. target_annotation_source: if the above target alignments were performed with bowtie2 or blast



## Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.tsv`.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
