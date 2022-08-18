
process COUNT_ANCHORS {

    tag "${fastq.simpleName}"
    label 'process_low'

    input:
    tuple val(fastq_id), path(fastq)

    output:
    path outfile, emit: seqs

    script:
    outfile = "counted_${fastq.simpleName}.txt"
    """
    zcat ${fastq} \\
        | sort -T ${workDir} \\
        | uniq -c > ${outfile}
    """
}
