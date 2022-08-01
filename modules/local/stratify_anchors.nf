
process STRATIFY_ANCHORS {

    label 'process_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    path counts
    val stratify_level
    val is_RNAseq

    output:
    path("stratified_*"), emit: seqs

    script:
    def is_RNAseq       = (is_RNAseq == true) ? "--is_RNAseq" : ""
    """
    for file in counted*txt
    do
        stratify_anchors.py \\
            --infile \${file} \\
            --stratify_level ${stratify_level}
    done
    """
}
