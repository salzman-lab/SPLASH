
process MAKE_FASTA {

    label 'process_low'

    input:
    path anchors

    output:
    path "*.fasta"       , emit: fasta

    script:
    """
    header=\$(cut -f1 ${anchors} | head -1)
    if [[ \${header} == "anchor" ]]
    then
        cut -f1 ${anchors} \
            | tail -n +2 \
            | sort \
            | uniq \
            | awk '{print ">"\$0"\\n"\$0}' \
            > anchor.fasta
    else
        cut -f1 ${anchors} \
            | sort \
            | uniq \
            | awk '{print ">"\$0"\\n"\$0}' \
            > anchor.fasta
    fi
    """
}
