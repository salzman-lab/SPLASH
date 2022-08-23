
process GET_ANCHOR_FASTA {

    tag "${samplesheet_id}"
    label 'process_low'

    input:
    tuple val(samplesheet_id), path(samplesheet), path(anchors_pvals)

    output:
    tuple val(samplesheet_id), path(samplesheet), path(outfile), emit: fasta

    script:
    outfile                  = "anchor_${samplesheet_id}.fasta"
    """
    cut -f1 ${anchors_pvals} \
        | tail -n +2 \
        | sort \
        | uniq \
        | awk '{print ">"\$0"\\n"\$0}' \
        > ${outfile}
    """
}
