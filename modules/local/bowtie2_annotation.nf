
process BOWTIE2_ANNOTATION {

    tag "${index_name}"
    label 'process_low'

    input:
    path anchor_fasta
    val index

    output:
    path "anchor_hits*tsv", emit: anchor_hits

    script:
    index_name = file(index).baseName
    index_hits = "anchor_hits_${index_name}.tsv"
    """
    rm -rf ${index_hits}
    echo -e "anchor\tanchor_hits_${index_name}\tanchor_hits_pos_${index_name}" >> ${index_hits}

    bowtie2 -f -x ${index} -U ${anchor_fasta} --quiet \\
        | sed '/^@/d' \\
        | cut -f1,3,4 \\
        | sort \\
        | uniq \\
        >> ${index_hits}
    """
}
