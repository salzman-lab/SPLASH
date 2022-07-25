
process DEDUP_CBC_UMI {

    label 'process_low'

    input:
    tuple val(fastq_id), path(fastq), val(group_id)

    output:
    tuple val(fastq_id), path(outfile), val(group_id), emit: seqs

    script:
    outfile         = "dedup_cbc_umi_${fastq_id}.txt.gz"
    """
    zcat ${fastq} \\
        | sort -u -k1,1 \\
        > ${outfile}
    """
}
