
process FETCH_ANCHORS {

    tag "${fastq_id}"
    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    tuple val(fastq_id), path(fastq)
    val is_10X
    val num_reads_first_pass
    val kmer_size
    val lookahead
    val anchor_mode
    val window_slide

    output:
    tuple val(fastq_id), path(outfile), emit: seqs

    script:
    def is_10X   = (is_10X == true) ? "--is_10X" : ""
    outfile      = "sequences_${fastq_id}.txt.gz"
    """
    fetch_anchors.py \\
        --infile ${fastq} \\
        --fastq_id ${fastq_id} \\
        --num_lines ${num_reads_first_pass} \\
        --kmer_size ${kmer_size} \\
        --lookahead ${lookahead} \\
        --anchor_mode ${anchor_mode} \\
        --window_slide ${window_slide} \\
        --outfile ${outfile} \\
        ${is_10X}
    """
}
