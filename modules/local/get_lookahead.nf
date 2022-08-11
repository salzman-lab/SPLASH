process GET_LOOKAHEAD {

    input:
    path fastq
    path samplesheet
    val kmer_size

    output:
    env lookahead, emit: lookahead

    script:
    """
    read_length=\$(zcat ${fastq} | head -n 2 | tail -n 1 | awk '{print length}' || true)

    d=\$(( (\${read_length} - 2 * ${kmer_size}) / 2 ))

    lookahead=\$(( \${d} > 0 ? \${d} : 0 ))
    """

}
