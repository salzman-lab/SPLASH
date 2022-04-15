
process NORM_SCORES {

    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 numpy=1.22.3" : null)


    input:
    path anchor_sample_scores
    path anchor_target_counts
    path samplesheet
    val use_std
    val kmer_size

    output:
    path ("*tsv")   , emit: norm_scores

    script:
    infile_hamming = "anchor_sample_scores_hamming.tsv"
    outfile_hamming = "anchor_norm_scores_hamming.tsv"
    infile_jaccard = "anchor_sample_scores_jaccard.tsv"
    outfile_jaccard = "anchor_norm_scores_jaccard.tsv"
    """
    norm_scores.py \\
        --anchor_sample_scores ${infile_hamming} \\
        --anchor_target_counts ${anchor_target_counts} \\
        --samplesheet ${samplesheet} \\
        --use_std ${use_std} \\
        --kmer_size ${kmer_size} \\
        --outfile_norm_scores ${outfile_hamming}

    norm_scores.py \\
        --anchor_sample_scores ${infile_jaccard} \\
        --anchor_target_counts ${anchor_target_counts} \\
        --samplesheet ${samplesheet} \\
        --use_std ${use_std} \\
        --kmer_size ${kmer_size} \\
        --outfile_norm_scores ${outfile_jaccard}
    """
}
