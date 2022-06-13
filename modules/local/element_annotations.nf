
process ELEMENT_ANNOTATIONS {

    label 'process_high_memory'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1" : null)

    input:
    path hits

    output:
    path outfile_ann_anchors, emit: annotated_anchors
    path outfile_ann_targets, emit: annotated_targets

    script:
    outfile_ann_anchors     = "annotated_anchors.tsv"
    outfile_ann_targets     = "annotated_targets.tsv"
    """
    element_annotations.py \\
        --outfile_ann_anchors ${outfile_ann_anchors} \\
        --outfile_ann_targets ${outfile_ann_targets}
    """
}
