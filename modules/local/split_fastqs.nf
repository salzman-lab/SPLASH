process SPLIT_FASTQS {

    tag "${fastq_id}"
    label 'process_medium'

    input:
    tuple val(fastq_id), path(fastq)
    val num_lines
    val num_chunk_lines

    output:
    path "*fastq.gz", emit: fastq

    script:
    """
    zcat ${fastq} \\
        | head -n ${num_lines} \\
        | split -l ${num_chunk_lines} --numeric 1 --numeric-suffixes - ${fastq_id}_ --additional-suffix .fastq \\
        || true

    for file in *fastq
    do
        gzip \${file}
    done
    """
}
