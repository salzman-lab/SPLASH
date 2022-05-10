
process MERGE_ANCHOR_SCORES {

    label 'process_low'

    input:
    path anchors

    output:
    path outfile    , emit: seqs

    script:
    outfile         = "anchors.tsv"
    """
    sort -k2n ${anchors} | head -n 5000 - | awk '{print \$1}' > ${outfile}
    """
}
