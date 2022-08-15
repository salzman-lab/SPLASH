
process COUNT_ANCHORS {

    tag "${fastq_id}"
    label 'process_low'

    input:
    tuple val(fastq_id), path(fastq)

    output:
    path outfile, emit: seqs

    script:
    outfile = "counted_${fastq_id}.txt"
    """
    zcat ${fastq} \\
        | sort -T ${workDir} \\
        | uniq -c > ${outfile}
    """
}
