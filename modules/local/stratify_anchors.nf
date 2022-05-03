
process STRATIFY_ANCHORS {

    label 'process_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    path counts
    val stratifiy_level

    output:
    path("stratified_*"), emit: seqs

    script:
    """
    for file in counted*txt
    do
        stratify_anchors.py \\
            --infile \${file} \\
            --stratifiy_level ${stratifiy_level}
    done
    """
}
