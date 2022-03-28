
process ADD_ANNOTATIONS {

    label 'process_medium'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    path anchor_hits_samplesheet
    path target_hits_samplesheet
    path anchor_scores

    output:
    path outfile_ann_anchor_scores
    path outfile_ann_targets

    script:
    outfile_ann_anchor_scores   = "anchor_scores.tsv"
    outfile_ann_targets         = "annotated_targets.tsv"
    """
    add_annotations.py \\
        --anchor_hits_samplesheet ${anchor_hits_samplesheet} \\
        --target_hits_samplesheet ${target_hits_samplesheet} \\
        --anchor_scores ${anchor_scores} \\
        --outfile_ann_anchor_scores ${outfile_ann_anchor_scores} \\
        --outfile_ann_targets ${outfile_ann_targets}
    """
}
