
process GET_ANCHORS_AND_SCORES {

    label 'process_medium'

    input:
    path counts
    path samplesheet
    val distance_type
    val max_targets
    val max_dist
    val bonfer
    val pval_threshold
    val run_type

    output:
    path outfile_scores     , emit: scores
    path outfile_anchors    , emit: anchors
    path "cmx*"             , emit: cmx

    script:
    file_id                 = counts.baseName
    outfile_scores          = "scores_${file_id}.tsv"
    outfile_anchors         = "anchors_${file_id}.tsv"
    """
    get_anchors.R \\
        ${counts} \\
        ${samplesheet} \\
        ${distance_type} \\
        ${max_targets} \\
        ${max_dist} \\
        ${bonfer} \\
        ${pval_threshold} \\
        ${run_type} \\
        ${outfile_scores} \\
        ${outfile_anchors}
    """
}
