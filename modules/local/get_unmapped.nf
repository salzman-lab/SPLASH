process GET_UNMAPPED {

    tag "${fastq_id}"
    label 'process_medium'

    input:
    tuple val(fastq_id), path(fastq), val(group_id)
    val index_bowtie

    output:
    tuple val(fastq_id), path("*unmapped.fastq.gz"), val(group_id) , emit: fastq


    script:
    """
    bowtie -n 0 --un ${fastq_id}_unmapped.fastq -x ${index_bowtie} ${fastq} > ${fastq_id}.sam
    gzip ${fastq_id}_unmapped.fastq
    rm ${fastq_id}.sam
    """
}
