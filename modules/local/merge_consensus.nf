
process MERGE_CONSENSUS {

    label 'process_low'

    input:
    path fastas

    output:
    path fasta  , emit: fasta

    script:
    fasta       = "merged_consensus.fasta"
    """
    cat *fasta > ${fasta}
    """
}
