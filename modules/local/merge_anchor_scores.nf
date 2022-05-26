
process MERGE_ANCHOR_SCORES {

    label 'process_low'

    input:
    path samplesheet
    val pval_threshold

    output:
    path outfile_scores     , emit: scores

    script:
    outfile_scores          = "anchors_pvals.tsv"
    """
    aggregate_pvals.py \\
        --infile ${samplesheet} \\
        --pval_threshold ${pval_threshold} \\
        --outfile_scores ${outfile_scores}
    """
}
