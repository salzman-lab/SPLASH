process GET_UNMAPPED {

    tag "${fastq_id}"
    label 'process_medium'

    input:
    tuple val(fastq_id), path(fastq)
    val index


    output:
    path "*fastq.gz", emit: fastq

    script:
    """
    for file in *fastq
    do
        bowtie2 -x ${index} -U ${fastq} --un-gz ${fastq_id}.fastq.gz > ${fastq_id}.sam
        rm ${fastq_id}.sam
    done
    """
}
