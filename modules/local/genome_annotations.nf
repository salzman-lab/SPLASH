
process GENOME_ANNOTATIONS {

    tag "${index_name}"
    label 'process_low'

    input:
    tuple path(fasta), path(genome_bam), path(transcriptome_bam)
    path gene_bed
    path exon_starts_bed
    path exon_ends_bed

    output:
    path "*tsv"

    script:
    fasta_name          = fasta.baseName
    anchor_hits         = "genome_annotations_anchor.tsv"
    target_hits         = "genome_annotations_target.tsv"
    """
    samtools view ${transcriptome_bam} \\
        | cut -f1,3 \\
        > \${fasta_name}_transcriptome_hit.txt

    samtools view ${transcriptome_bam} \\
        | cut -f1,5 \\
        > \${fasta_name}_transcriptome_mapq.txt

    bedtools bamtobed -i ${genome_bam} \\
        | sort -k1,1 -k2,2n \\
        > \${fasta_name}_genome.bed

    if [[ \$(wc -l \${fasta_name}_genome.bed | awk '{print \$1}') -gt 0 ]]
    then
        cut -f4,6 \${fasta_name}_genome.bed \\
            | sort \\
            > \${fasta_name}_genome_strand.txt
        cut -f4,5 \${fasta_name}_genome.bed \\
            | sort \\
            > \${fasta_name}_genome_mapq.txt
        bedtools intersect -a \${fasta_name}_genome.bed -b ${gene_bed} -wb \\
            | cut -f1-4,10,6 \\
            | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
            | cut -f4,6 \\
            | sort \\
            > \${fasta_name}_genome_genes.txt
        bedtools closest -a \${fasta_name}_genome.bed -b ${exon_starts_bed} -D ref -id -t first \\
            | cut -f4,13 \\
            | sort \\
            > \${fasta_name}_genome_upstream_exon_starts.txt
        bedtools closest -a \${fasta_name}_genome.bed -b ${exon_ends_bed} -D ref -id -t first \\
            | cut -f4,13 \\
            | sort \\
            > \${fasta_name}_genome_upstream_exon_ends.txt
        bedtools closest -a \${fasta_name}_genome.bed -b ${exon_starts_bed} -D ref -iu -t first \\
            | cut -f4,13 \\
            | sort \\
            > \${fasta_name}_genome_downstream_exon_starts.txt
        bedtools closest -a \${fasta_name}_genome.bed -b ${exon_ends_bed} -D ref -iu -t first \\
            | cut -f4,13 \\
            | sort \\
            > \${fasta_name}_genome_downstream_exon_ends.txt
    fi

    merge_genome_annotations.py \\
        --fasta_name ${fasta_name}
    """
}
