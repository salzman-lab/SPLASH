
process COUNT_ANCHORS {

    label 'process_medium'

    input:
    path bam

    output:
    path "*tsv"

    script:
    bam_name    = bam.baseName
    outfile     = "ann_splices_${bam_name}.txt"
    """
    
    """
}
