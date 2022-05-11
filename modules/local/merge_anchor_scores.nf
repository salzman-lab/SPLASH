
process MERGE_ANCHOR_SCORES {

    label 'process_low'

    input:
    path all_anchors

    output:
    path anchors        , emit: anchors
    path anchors_pvals  , emit: anchors_pvals

    script:
    anchors             = "anchors.tsv"
    anchors_pvals       = "anchors_pvals.tsv"
    """
    sort -k2n ${all_anchors} \
        | head -n 5000 - \
        > ${anchors_pvals} || true

    awk '{print \$1}' ${anchors_pvals} \
        > ${anchors} || true
    """
}
