
process MERGE_ABUNDANT_ANCHORS {

    label 'process_medium'

    input:
    path counts
    val num_decoy_anchors

    output:
    path outfile    , emit: seqs

    script:
    outfile         = "most_abundant_anchors.tsv"
    """
    cat counts* \\
        | sort -k1nr \\
        | head -n ${num_decoy_anchors} \\
        | cut -f2 -d" " \\
        > ${outfile} \\
        || true
    """
}
