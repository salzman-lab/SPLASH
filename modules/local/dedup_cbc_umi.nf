process DEDUP_CBC_UMI {

    tag "${id}"
    label 'process_low'

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path(outfile), emit: seqs

    script:
    outfile         = "dedup_cbc_umi_${id}.txt.gz"
    """
    zcat ${fastq} \\
        | sort -u -k1,1 \\
        > dedup_cbc_umi_${id}.txt
        
    gzip dedup_cbc_umi_${id}.txt
    """
}
