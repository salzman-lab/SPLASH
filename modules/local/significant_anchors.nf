
process SIGNIFICANT_ANCHORS {

    label 'process_high'

    input:
    path anchors
    val fdr_threshold
    val all_anchors_pvals_file

    output:
    path outfile_scores         , emit: scores
    path "all*tsv"              , emit: all_scores

    script:
    def all_anchors_pvals_file  = all_anchors_pvals_file   ? "--all_anchors_pvals_file" : ""

    outfile_scores              = "anchors_pvals.tsv"
    """
    significant_anchors.py \\
        --fdr_threshold ${fdr_threshold} \\
        --outfile_scores ${outfile_scores} \\
        ${all_anchors_pvals_file}
    """
}
