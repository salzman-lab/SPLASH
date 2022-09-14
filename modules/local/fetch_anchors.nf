
process FETCH_ANCHORS {

    tag "${id}"
    label 'process_very_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    tuple val(id), path(fastq)
    val is_10X
    val num_fastq_reads
    val kmer_size
    val lookahead
    val anchor_mode
    val window_slide

    output:
    tuple val(id), path(outfile), emit: seqs

    script:
    def is_10X  = (is_10X == true) ? "--is_10X" : ""
    fastq_id    = "${id}_${fastq.simpleName}"
    outfile     = "sequences_${fastq_id}.txt.gz"
    """
    fetch_anchors.py \\
        --infile ${fastq} \\
        --id ${id} \\
        --num_fastq_reads ${num_fastq_reads} \\
        --kmer_size ${kmer_size} \\
        --lookahead ${lookahead} \\
        --anchor_mode ${anchor_mode} \\
        --window_slide ${window_slide} \\
        --outfile ${outfile} \\
        ${is_10X}
    """
}
