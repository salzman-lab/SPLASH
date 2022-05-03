
process GENOME_ANNOTATIONS {

    tag "${index_name}"
    label 'process_low'

    input:
    path anchor_bam
    path target_bam
    path gene_bed
    path exon_starts_bed
    path exon_ends_bed

    output:
    path anchor_hits    , emit: anchor_annotations, optional: true
    path target_hits    , emit: target_annotations, optional: true

    script:
    anchor_hits         = "genome_annotations_anchors.tsv"
    target_hits         = "genome_annotations_targets.tsv"
    """
    bedtools bamtobed -i ${anchor_bam} | sort -k1,1 -k2,2n > anchors_genome.bed
    if [[ \$(wc -l anchors_genome.bed) -gt 0 ]]
    then
        cut -f4,6 anchors_genome.bed | sort > anchors_strand.txt
        bedtools intersect -a anchors_genome.bed -b ${gene_bed} | cut -f1-4,10,6 | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse | cut -f4,6 | sort > anchors_genes.txt
        bedtools closest -a anchors_genome.bed -b ${exon_starts_bed} -D ref -t first | cut -f4,13 | sort > anchors_exon_starts_distances.txt
        bedtools closest -a anchors_genome.bed -b ${exon_ends_bed} -D ref -t first | cut -f4,13 | sort > anchors_exon_ends_distances.txt
        join anchors_strand.txt anchors_genes.txt | join - anchors_exon_starts_distances.txt | join - anchors_exon_ends_distances.txt | sed 's/ /\t/g' > ${anchor_hits}
    fi

    bedtools bamtobed -i ${target_bam} | sort -k1,1 -k2,2n > targets_genome.bed
    if [[ \$(wc -l targets_genome.bed) -gt 0 ]]
    then
        cut -f4,6 targets_genome.bed | sort > targets_strand.txt
        bedtools intersect -a targets_genome.bed -b ${gene_bed} | cut -f1-4,10,6 | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse | cut -f4,6 | sort > targets_genes.txt
        bedtools closest -a targets_genome.bed -b ${exon_starts_bed} -D ref -t first | cut -f4,13 | sort > targets_exon_starts_distances.txt
        bedtools closest -a targets_genome.bed -b ${exon_ends_bed} -D ref -t first | cut -f4,13 | sort > targets_exon_ends_distances.txt
        join targets_strand.txt targets_genes.txt | join - targets_exon_starts_distances.txt | join - targets_exon_ends_distances.txt | sed 's/ /\t/g' > ${target_hits}
    fi
    """
}
