
process SPLICING_ANNOTATIONS {

    label 'process_medium'
    conda (params.enable_conda ? "bioconda::samtools bioconda::bedtools=2.30.0 conda-forge::python=3.9.5 pandas=1.4.3 numpy=1.22.3" : null)

    input:
    path bam
    path unmapped_fasta
    path gene_bed
    path fasta
    path genome_annotations_anchors

    output:
    path outfile_annotations    , emit: consensus_called_exons  , optional: true
    path "*tsv"                 , emit: tsv                     , optional: true
    path "called_exons.bed"     , emit: bed                     , optional: true
    path fasta                  , emit: fasta                   , optional: true
    path "consensus_genes.txt"  , emit: consenus_genes          , optional: true
    path "*bam*"                , emit: bam_bai                 , optional: true

    script:
    outfile_unmapped            = "unmapped_consensus_sequences.tsv"
    outfile_annotations         = "consensus_called_exons.tsv"
    """
    ## sort bam and index for later
    samtools sort ${bam} > sorted_${bam}
    samtools index sorted_${bam}

    ## get reported alignments
    samtools view ${bam} \\
        | cut -f1,10 \\
        >> reported_alignments.txt

    ## get called exons
    bedtools bamtobed -split -i ${bam} \\
        | sed '/^chr/!d' \\
        | sort -k1,1 -k2,2n \\
        > called_exons.bed

    ## get called exons start and end positions
    awk -v OFS='\\t' '{print \$1,\$2-1,\$2+1,\$4,"called_exon_start",\$6"\\n"\$1,\$3-1,\$3+1,\$4,"called_exon_end",\$6}' called_exons.bed \\
        | sort -k1,1 -k2,2n \\
        | awk -v OFS='\\t' '{if (\$2 < 0) print \$1,0,\$3,\$4,\$5,\$6; else print \$1,\$2,\$3,\$4,\$5,\$6}' \\
        > positions_called_exons.bed

    ## get consensus genes
    bedtools intersect -a called_exons.bed -b ${gene_bed} -wb -loj \\
        | cut -f 4,10 \\
        | bedtools groupby -g 1 -c 2 -o distinct \\
        > consensus_genes.txt

    ## add consensus and anchor gene columns
    splicing_annotations.py \\
        --unmapped_fasta ${unmapped_fasta} \\
        --fasta ${fasta} \\
        --ann_called_exons positions_called_exons.bed \\
        --genome_annotations_anchors ${genome_annotations_anchors} \\
        --consensus_genes consensus_genes.txt \\
        --reported_alignments reported_alignments.txt \\
        --outfile_unmapped ${outfile_unmapped} \\
        --outfile_annotations ${outfile_annotations}
    """
}
