
process ADDITIONAL_SUMMARY {

    label 'process_low'

    input:
    path consensus_genes
    path summary

    output:
    path outfile                , emit: tsv

    script:
    outfile                     = "additional_summary.tsv"
    """
    additional_summary.py \\
        --summary ${summary} \\
        --consensus_genes ${consensus_genes} \\
        --outfile ${outfile}
    """
}
