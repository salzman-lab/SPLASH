
process ADDITIONAL_SUMMARY {

    label 'process_high_memory'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.3" : null)

    input:
    path consensus_called_exons
    path summary

    output:
    path outfile                , emit: tsv

    script:
    outfile                     = "additional_summary.tsv"
    """
    additional_summary.py \\
        --summary ${summary} \\
        --consensus_called_exons ${consensus_called_exons} \\
        --outfile ${outfile}
    """
}
