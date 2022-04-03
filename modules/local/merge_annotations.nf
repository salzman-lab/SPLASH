
process MERGE_ANNOTATIONS {

    label 'process_medium'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1" : null)

    input:
    path anchor_hits_samplesheet
    path target_hits_samplesheet

    output:
    path outfile_ann_anchors
    path outfile_ann_targets

    script:
    outfile_ann_anchors     = "annotated_anchors.tsv"
    outfile_ann_targets     = "annotated_targets.tsv"
    """
    merge_annotations.py \\
        --anchor_hits_samplesheet ${anchor_hits_samplesheet} \\
        --target_hits_samplesheet ${target_hits_samplesheet} \\
        --outfile_ann_anchors ${outfile_ann_anchors} \\
        --outfile_ann_targets ${outfile_ann_targets}
    """
}
