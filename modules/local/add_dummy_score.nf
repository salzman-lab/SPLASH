
process ADD_DUMMY_SCORE {

    label 'process_low'

    input:
    path anchors

    output:
    path outfile    , emit: anchors_pvals

    script:
    outfile         = "anchors_pvals.tsv"
    """
    awk -v OFS='\\t' '{print \$1, "NaN"}' ${anchors} \\
        > ${outfile}
    """
}
