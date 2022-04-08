process GET_UNMAPPED {

    tag "${fastq_id}"
    label 'process_medium'

    input:
    tuple val(fastq_id), path(fastq)
    val index_bowtie

    output:
    tuple val(fastq_id), path("*unmapped.fastq.gz")        , emit: fastq


    script:
    """
    bowtie2 -x ${index_bowtie} -U ${fastq} --un-gz ${fastq_id}_unmapped.fastq.gz > ${fastq_id}.sam
    rm ${fastq_id}.sam
    """
}
