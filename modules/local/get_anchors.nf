
process GET_ANCHORS {

    input:
    path ch_split_fastqs
    val chunk_size
    path samplesheet
    val target_threshold
    val anchor_counts_threshold
    val anchor_freeze_threshold
    val anchor_mode
    val window_slide
    val read_len

    output:
    path "anchors.tsv", emit: anchors

    script:
    outfile = "anchors.tsv"
    """
    get_anchors.py \\
        --n_iterations 100 \\
        --max_reads ${chunk_size} \\
        --samplesheet ${samplesheet} \\
        --target_threshold ${} \\
        --anchor_counts_threshold ${} \\
        --anchor_freeze_threshold ${} \\
        --anchor_mode ${} \\
        --window_slide ${} \\
        --read_len ${} \\
        --outfile ${outfile}
    """
}
