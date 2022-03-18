process SPLIT_FASTQS {

    input:
    tuple val(fastq_id), path(fastq)
    val num_lines
    val num_chunk_lines

    output:
    path "*fastq.gz", emit: fastq

    script:
    """
    zcat ${fastq} \\
        | head -${num_lines} \\
        | split -l ${num_chunk_lines} --numeric-suffixes - ${fastq_id}_ --additional-suffix .fastq

    for file in *fastq
    do
        gzip \${file}
    done
    """
}
