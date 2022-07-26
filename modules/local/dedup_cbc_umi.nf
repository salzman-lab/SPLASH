
process DEDUP_CBC_UMI {

    label 'process_low'

    input:
    tuple val(fastq_id), path(fastq)

    output:
    tuple val(fastq_id), path(outfile), emit: seqs

    script:
    outfile         = "dedup_cbc_umi_${fastq_id}.txt.gz"
    """
    zcat ${fastq} \\
        | sort -u -k1,1 \\
        > dedup_cbc_umi_${fastq_id}.txt

    gzip dedup_cbc_umi_${fastq_id}.txt
    """
}
