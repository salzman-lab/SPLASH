
process MERGE_ABUNDANT_ANCHORS {

    label 'process_medium'

    input:
    path counts

    output:
    path outfile    , emit: seqs

    script:
    outfile         = "most_abundant_anchors.tsv"
    """
    cat counts* \\
        | sort -k1nr \\
        | head -n 150 \\
        | cut -f2 -d" " \\
        > ${outfile} \\
        || true
    """
}
