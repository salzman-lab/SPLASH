
process ADDITIONAL_SUMMARY {

    label 'process_high_memory'

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
