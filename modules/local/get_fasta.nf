
process GET_FASTA {

    label 'process_low'

    input:
    path anchor_target_counts

    output:
    path "anchors.fasta"  , emit: anchor_fasta
    path "targets.fasta"  , emit: target_fasta

    script:
    """
    cut -f1 ${anchor_target_counts} \
        | tail -n +2 \
        | sort \
        | uniq \
        | awk '{print ">"\$0"\\n"\$0}' > anchors.fasta

    cut -f2 ${anchor_target_counts} \
        | tail -n +2 \
        | sort \
        | uniq \
        | awk '{print ">"\$0"\\n"\$0}' > targets.fasta
    """
}
