
process MERGE_ANCHOR_HITS {

    label 'process_low'

    input:
    path anchor_hits_samplesheet
    path anchor_scores

    output:
    path outfile

    script:
    outfile = "anchor_scores.tsv"
    """
    merge_anchor_hits.py \\
        --samplesheet ${anchor_hits_samplesheet} \\
        --summary_table ${anchor_scores} \\
        --outfile ${outfile}
    """
}
