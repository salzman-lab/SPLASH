
process COUNT_ANCHORS {

    tag "${id}"
    label 'process_low'

    input:
    tuple val(id), path(fastq)

    output:
    path outfile, emit: seqs

    script:
    outfile = "counted_${fastq.simpleName}.txt"
    """
    zcat ${fastq} \\
        | sort \\
        | uniq -c > ${outfile}
    """
}
