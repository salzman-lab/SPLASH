
process MERGE_ABUNDANT_ANCHORS {

    label 'process_medium'

    input:
    path counts

    output:
    path outfile    , emit: seqs

    script:
    outfile         = "merged_abundant_anchor_counts.tsv.gz"
    """
    cat abundant*gz > ${outfile}
    """
}
