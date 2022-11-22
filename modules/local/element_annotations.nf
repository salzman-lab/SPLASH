
process ELEMENT_ANNOTATIONS {

    label 'process_medium'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.3" : null)

    input:
    path hits

    output:
    path outfile_ann_anchors, emit: anchors

    script:
    outfile_ann_anchors     = "element_annotations_anchors.tsv"
    """
    element_annotations.py \\
        --outfile_ann_anchors ${outfile_ann_anchors}
    """
}
