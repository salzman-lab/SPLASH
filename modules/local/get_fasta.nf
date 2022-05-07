
process GET_FASTA {

    label 'process_low'

    input:
    path anchor_target_counts

    output:
    path "anchor.fasta"  , emit: anchor_fasta
    path "target.fasta"  , emit: target_fasta
    path "anchor.txt"    , emit: anchors
    path "target.txt"    , emit: targets

    script:
    """
    cut -f1 ${anchor_target_counts} \
        | tail -n +2 \
        | sort \
        | uniq \
        > anchor.txt
    awk '{print ">"\$0"\\n"\$0}' anchor.txt \
        > anchor.fasta

    cut -f2 ${anchor_target_counts} \
        | tail -n +2 \
        | sort \
        | uniq \
        > target.txt

    awk '{print ">"\$0"\\n"\$0}' target.txt \
        > target.fasta
    """
}
