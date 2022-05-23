process GET_READ_LENGTH {

    input:
    tuple val(fastq_id), path(fastq), val(group_id)

    output:
    env read_length, emit: read_length

    script:
    """
    read_length=\$(zcat ${fastq} | head -n 2 | tail -n 1 | awk '{print length}' || true)
    """

}
