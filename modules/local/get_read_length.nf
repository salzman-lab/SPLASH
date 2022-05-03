process GET_READ_LENGTH {

    input:
    tuple val(fastq_id), path(fastq), val(group_id)

    output:
    env read_length, emit: read_length

    script:
    if (fastq.getExtension() == 'gz')
        """
        zcat ${fastq} > tmp.fastq
        read_length=\$(head -n 2 tmp.fastq | tail -n 1 | awk '{print length}')
        rm tmp.fastq
        """
    else
        """
        read_length=\$(head -n 2 ${fastq} | tail -n 1 | awk '{print length}')
        """
}
