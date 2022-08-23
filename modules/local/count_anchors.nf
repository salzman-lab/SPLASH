
process COUNT_ANCHORS {

    tag "${id}"
    label 'process_low'

    input:
    tuple val(id), path(fastq)

    output:
    path outfile, emit: seqs

    script:
    fastq_id    = "${id}_${fastq.simpleName}"
    outfile     = "counted_${fastq_id}.txt"
    """
    zcat ${fastq} \\
        | sort -T ${workDir} \\
        | uniq -c > ${outfile}
    """
}
