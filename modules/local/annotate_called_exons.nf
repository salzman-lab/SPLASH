
process ANNOTATE_CALLED_EXONS {

    label 'process_low'
    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)

    input:
    path bam
    path ann_AS_gtf
    path fasta
    path genome_annotations_anchors

    output:
    path "*tsv"                 , emit: tsv
    path "called_exons.bed"     , emit: bed
    path fasta                  , emit: fasta

    script:
    outfile                     = "consensus_called_exons.tsv"
    """
    ## get called exons
    bedtools bamtobed -split -i ${bam} \\
        > called_exons.bed

    ## get called exons start and end positions
    awk -v OFS='\t' '{print \$1,\$2-1,\$2+1,\$4,"called_exon_start",\$6"\\n"\$1,\$3-1,\$3+1,\$4,"called_exon_end",\$6}' called_exons.bed \\
        | sort -k1,1 -k2,2n \\
        > positions_called_exons.bed

    ## annotate called exon start and ends
    bedtools intersect -a positions_called_exons.bed -b ${ann_AS_gtf} -loj -wb -sorted \\
        | cut -f1-5,10-13 \\
        | bedtools groupby -g 1,2,3,4,5, -c 6,7,8,9 -o distinct,distinct,max,max \\
        > annotated_positions_called_exons.bed

    ## add consensus and anchor gene columns
    merge_consensus_splice_annotation.py \\
        --infile annotated_positions_called_exons.bed \\
        --fasta ${fasta} \\
        --genome_annotations_anchors ${genome_annotations_anchors} \\
        --outfile ${outfile}
    """
}
