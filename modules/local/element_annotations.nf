
process ELEMENT_ANNOTATIONS {

    label 'process_high_memory'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.3" : null)

    input:
    tuple val(samplesheet_id), path(hits)

    output:
    tuple val(samplesheet_id), path(outfile_ann_anchors), emit: annotated_anchors, optional: true

    script:
    outfile_ann_anchors         = "element_annotations_anchors_${samplesheet_id}.tsv"
    """
    element_annotations.py \\
        --outfile_ann_anchors ${outfile_ann_anchors}
    """
}
