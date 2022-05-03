
process GENOME_ANNOTATIONS {

    tag "${index_name}"
    label 'process_low'

    input:
    path anchor_genome_bam
    path target_genome_bam
    path anchor_trans_bam
    path target_trans_bam
    path anchors
    path targets
    path gene_bed
    path exon_starts_bed
    path exon_ends_bed

    output:
    path anchor_hits
    path target_hits

    script:
    anchor_hits         = "genome_annotations_anchors.tsv"
    target_hits         = "genome_annotations_targets.tsv"
    """
    ## Process for anchors
    bedtools bamtobed -i ${anchor_genome_bam} \\
        | sort -k1,1 -k2,2n \\
        > anchors_genome.bed
    if [[ \$(wc -l anchors_genome.bed | awk '{print \$1}') -gt 0 ]]
    then
        cut -f4,6 anchors_genome.bed \\
            | sort \\
            > anchors_strand.txt
        bedtools intersect -a anchors_genome.bed -b ${gene_bed} -wb \\
            | cut -f1-4,10,6 \\
            | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
            | cut -f4,6 \\
            | sort \\
            > anchors_genes.txt
        bedtools closest -a anchors_genome.bed -b ${exon_starts_bed} -D ref -t first \\
            | cut -f4,13 \\
            | sort \\
            > anchors_exon_starts_distances.txt
        bedtools closest -a anchors_genome.bed -b ${exon_ends_bed} -D ref -t first \\
            | cut -f4,13 \\
            | sort \\
            > anchors_exon_ends_distances.txt
        samtools view ${anchor_trans_bam} \\
            | cut -f1,3 \\
            > anchors_transcriptome.txt
        rm -rf ${anchor_hits}
        echo -e "anchor\tstrand\tgene\tdistance_exon_start\tdistance_exon_end\ttranscript" >> ${anchor_hits}
        join <(sort ${anchors}) anchors_strand.txt -a1 -a2 -e- \\
            | join - anchors_genes.txt -a1 -a2 -e- \\
            | join - anchors_exon_starts_distances.txt -a1 -a2 -e- \\
            | join - anchors_exon_ends_distances.txt -a1 -a2 -e- \\
            | join - anchors_transcriptome.txt -a1 -a2 -e- \\
            | sed 's/ /\t/g' \\
            >> ${anchor_hits}
    fi

    ## Process for targets
    bedtools bamtobed -i ${target_genome_bam} \\
        | sort -k1,1 -k2,2n \\
        > targets_genome.bed
    if [[ \$(wc -l targets_genome.bed | awk '{print \$1}') -gt 0 ]]
    then
        cut -f4,6 targets_genome.bed \\
            | sort \\
            > targets_strand.txt
        bedtools intersect -a targets_genome.bed -b ${gene_bed} -wb \\
            | cut -f1-4,10,6 \\
            | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
            | cut -f4,6 \\
            | sort \\
            > targets_genes.txt
        bedtools closest -a targets_genome.bed -b ${exon_starts_bed} -D ref -t first \\
            | cut -f4,13 \\
            | sort \\
            > targets_exon_starts_distances.txt
        bedtools closest -a targets_genome.bed -b ${exon_ends_bed} -D ref -t first \\
            | cut -f4,13 \\
            | sort \\
            > targets_exon_ends_distances.txt
        samtools view ${target_trans_bam} \\
            | cut -f1,3 \\
            > targets_transcriptome.txt
        rm -rf ${target_hits}
        echo -e "target\tstrand\tgene\tdistance_exon_start\tdistance_exon_end\ttranscript" >> ${target_hits}
        join <(sort ${targets}) targets_strand.txt -a1 -a2 -e- \\
            | join - targets_genes.txt -a1 -a2 -e- \\
            | join - targets_exon_starts_distances.txt -a1 -a2 -e- \\
            | join - targets_exon_ends_distances.txt -a1 -a2 -e- \\
            | join - targets_transcriptome.txt -a1 -a2 -e- \\
            | sed 's/ /\t/g' \\
            >> ${target_hits}
    fi
    """
}
