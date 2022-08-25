
process GENOME_ANNOTATIONS {

    tag "${samplesheet_id}"
    publishDir (
        path: {"${params.outdir}/${samplesheet_id}/genome_annotations"},
        mode: "copy",
        pattern: "*.tsv")
    label 'process_high'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6 pandas' : null)

    input:
    tuple val(samplesheet_id), path(samplesheet), path(fasta), path(genome_bam), path(transcriptome_bam)
    path gene_bed

    output:
    tuple val(samplesheet_id), path(anchor_hits) , emit: annotated_anchors, optional: true

    script:
    anchor_hits = "genome_annotations_anchors.tsv"
    """
    samtools view ${transcriptome_bam} \\
        | cut -f1,3 \\
        > ${samplesheet_id}_transcriptome_hit.txt

    samtools view ${transcriptome_bam} \\
        | cut -f1,5 \\
        >  ${samplesheet_id}_transcriptome_mapq.txt

    bedtools bamtobed -i ${genome_bam} \\
        | sort -k1,1 -k2,2n \\
        >  ${samplesheet_id}_genome.bed

    if [[ \$(wc -l  ${samplesheet_id}_genome.bed | awk '{print \$1}') -gt 0 ]]
    then
        cut -f4,6  ${samplesheet_id}_genome.bed \\
            | sort \\
            >  ${samplesheet_id}_genome_strand.txt
        cut -f4,5 ${samplesheet_id}_genome.bed \\
            | sort \\
            >  ${samplesheet_id}_genome_mapq.txt
        bedtools intersect -a  ${samplesheet_id}_genome.bed -b ${gene_bed} -wb \\
            > intersect_genome_genomes.bed

        if [[ \$(wc -l  intersect_genome_genomes.bed | awk '{print \$1}') -gt 0 ]]
        then
            cut -f1-4,10,6 intersect_genome_genomes.bed \\
                | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
                | cut -f4,6 \\
                | sort \\
                >  ${samplesheet_id}_genome_genes.txt
        fi
    fi

    genome_annotations.py \\
        --fasta_name ${samplesheet_id}
    """
}
