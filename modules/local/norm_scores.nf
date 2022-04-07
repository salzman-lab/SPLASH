
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
    path outfile_norm_scores   , emit: norm_scores

    script:
    outfile_norm_scores        = "anchor_norm_scores.tsv"
    """
    norm_scores.py \\
        --anchor_sample_scores ${anchor_sample_scores} \\
        --anchor_target_counts ${anchor_target_counts} \\
        --samplesheet ${samplesheet} \\
        --use_std ${use_std} \\
        --kmer_size ${kmer_size} \\
        --outfile_norm_scores ${outfile_norm_scores}
    """
}
