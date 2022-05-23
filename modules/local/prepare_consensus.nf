
process PREPARE_CONSENSUS {

    label 'process_low'

    input:
    path fasta

    output:
    path fasta  , emit: fasta

    script:
    fasta_id    = fasta.baseName
    """
    sed -i "s/>/>${fasta_id}____/g" ${fasta}
    """
}
