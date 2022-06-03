
process SIGNIFICANT_ANCHORS {

    label 'process_high'

    input:
    path anchors
    val fdr_threshold

    output:
    path outfile_scores     , emit: scores

    script:
    outfile_scores          = "anchors_pvals.tsv"
    """
    significant_anchors.py \\
        --fdr_threshold ${fdr_threshold} \\
        --outfile_scores ${outfile_scores}
    """
}
