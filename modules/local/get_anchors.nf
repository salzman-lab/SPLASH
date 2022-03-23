
process GET_ANCHORS {

    label 'process_high'

    input:
    path ch_split_fastqs
    val n_iterations
    val chunk_size
    val kmer_size
    path samplesheet
    val target_counts_threshold
    val anchor_counts_threshold
    val anchor_freeze_threshold
    val anchor_score_threshold
    val anchor_mode
    val c_type
    val window_slide
    val looklength
    val num_keep_anchors
    val use_std
    val compute_target_distance

    output:
    path "anchors.tsv", emit: anchors
    path "ignorelist*", emit: ignorelist
    path "*log"

    script:
    outfile = "anchors.tsv"
    """
    get_anchors.py \\
        --n_iterations ${n_iterations} \\
        --max_reads ${chunk_size} \\
        --kmer_size ${kmer_size} \\
        --samplesheet ${samplesheet} \\
        --target_counts_threshold ${target_counts_threshold} \\
        --anchor_counts_threshold ${anchor_counts_threshold} \\
        --anchor_freeze_threshold ${anchor_freeze_threshold} \\
        --anchor_score_threshold ${anchor_score_threshold} \\
        --anchor_mode ${anchor_mode} \\
        --window_slide ${window_slide} \\
        --c_type ${c_type} \\
        --looklength ${looklength} \\
        --num_keep_anchors ${num_keep_anchors} \\
        --use_std ${use_std} \\
        --compute_target_distance ${compute_target_distance} \\
        --outfile ${outfile}
    """
}
