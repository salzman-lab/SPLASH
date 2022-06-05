
process COUNT_ANCHORS {

    tag "${fastq_id}"
    label 'process_low'

    input:
    tuple val(fastq_id), path(fastq), val(group_id)

    output:
    path outfile, emit: seqs

    script:
    outfile = "counted_${fastq_id}.txt"
    """
    zcat ${fastq} \\
        | sort \\
        | uniq -c > ${outfile}
    """
}
