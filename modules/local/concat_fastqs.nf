process CONCAT_FASTQS {

    label 'process_low'

    input:
    tuple val(id), path(fastqs)

    output:
    tuple val(id), path(seqs) , emit: seqs

    script:
    seqs                      = "merged_${id}.txt.gz"
    """
    rm -rf ${seqs}
    cat extracted*txt.gz > ${seqs}
    """
}
