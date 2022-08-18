process CONCAT_FASTQS {

    label 'process_low'

    input:
    tuple val(id), path(fastqs)

    output:
    tuple val(id), path(fastq) , emit: fastqs

    script:
    fastq                           = "merged_${id}.fastq.gz"
    """
    rm -rf ${fastq}
    cat *fastq.gz > ${fastq}
    """
}
