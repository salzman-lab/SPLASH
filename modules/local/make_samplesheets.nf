
process MAKE_SAMPLESHEETS {

    tag "samplesheets"
    label "process_low"

    input:
    path base_samplesheet
    val is_10X
    val run_unsupervised_pvals

    output:
    path "*csv"         , emit: samplesheets    , optional: true

    script:
    def is_10X                  = (is_10X == true)                  ? "--is_10X" : ""
    def run_unsupervised_pvals  = (run_unsupervised_pvals == true)  ? "--run_unsupervised_pvals --outfile unsupervised_pvalues_samplesheet.csv" : ""
    """
    make_pairwise_samplesheets.py \\
        --samplesheet ${base_samplesheet} \\
        ${is_10X} \\
        ${run_unsupervised_pvals}
    """
}
