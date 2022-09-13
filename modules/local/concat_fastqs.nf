process CONCAT_FASTQS {

    label 'process_low'

    input:
    tuple val(id), path(fastqs)
    val num_fastq_reads

    output:
    tuple val(id), path(merged_fastq)   , emit: fastqs

    script:
    merged_fastq                        = "merged_${id}.fastq.gz"
    """
    ## Expand number of reads by 4 for FASTQ files
    num_reads=\$(( 4 * ${num_fastq_reads} ))

    ## Clear out file for addition by appending
    rm -rf ${merged_fastq}

    ## For each FASTQ file, sample the number of reads and append to the channel FASTQ file
    for fastq in *fastq.gz
    do
        zcat \${fastq} \\
            | head -n \${num_reads} - \\
            >> merged_${id}.fastq \\
            || true
    done

    ## gzip
    gzip merged_${id}.fastq

    """
}
