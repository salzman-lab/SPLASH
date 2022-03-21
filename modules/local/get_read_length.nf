process GET_READ_LENGTH {

    input:
    tuple val(fastq_id), path(fastq

    output:
    env read_length, emit: read_length

    script:
    """
    zcat ${fastq} > tmp.fastq
    read_length=\$(head -n 1 tmp.fastq | awk '{print length}')
    """
}
