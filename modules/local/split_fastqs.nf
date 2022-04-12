process SPLIT_FASTQS {

    tag "${fastq_id}"
    label 'process_medium'

    input:
    tuple val(fastq_id), path(fastq)
    val num_lines
    val num_chunk_lines
    val n_iterations

    output:
    path "*fastq.gz", emit: fastq

    script:
    suffix_length = n_iterations.toString().length()
    """
    zcat ${fastq} \\
        | head -n ${num_lines} \\
        | split -l ${num_chunk_lines} --numeric=1 --numeric-suffixes - ${fastq_id}_ --additional-suffix .fastq --suffix-length=${suffix_length} \\
        || true

    for file in *fastq
    do
        if [[ ! -f \${file}.gz ]]
        then
            gzip \${file}
        fi
    done
    """
}
