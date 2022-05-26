
process SIGNIFICANT_ANCHORS {

    label 'process_low'

    input:
    path samplesheet
    val pval_threshold

    output:
    path outfile_scores     , emit: scores

    script:
    outfile_scores          = "anchors_pvals.tsv"
    """
    significant_anchors.py \\
        --infile ${samplesheet} \\
        --pval_threshold ${pval_threshold} \\
        --outfile_scores ${outfile_scores}
    """
}
