
process GET_ANCHORS {

    label 'process_medium'

    input:
    path ch_split_fastqs
    val n_iterations
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
    path ".*log"

    script:
    outfile = "anchors.tsv"
    """
    get_anchors.py \\
        --n_iterations ${n_iterations} \\
        --max_reads ${chunk_size} \\
        --samplesheet ${samplesheet} \\
        --target_threshold ${target_threshold} \\
        --anchor_counts_threshold ${anchor_counts_threshold} \\
        --anchor_freeze_threshold ${anchor_freeze_threshold} \\
        --anchor_mode ${anchor_mode} \\
        --window_slide ${window_slide} \\
        --read_len ${read_len} \\
        --outfile ${outfile}
    """
}
