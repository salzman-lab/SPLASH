
process GENOME_ANNOTATIONS {

    tag "${fasta_name}"
    label 'process_medium'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6 pandas' : null)

    input:
    tuple path(fasta), path(genome_bam), path(transcriptome_bam)
    path gene_bed

    output:
    path anchor_hits    , emit: anchors

    script:
    fasta_name          = fasta.baseName
    anchor_hits         = "genome_annotations_anchor.tsv"
    """
    samtools view ${transcriptome_bam} \\
        | cut -f1,3 \\
        > ${fasta_name}_transcriptome_hit.txt

    samtools view ${transcriptome_bam} \\
        | cut -f1,5 \\
        >  ${fasta_name}_transcriptome_mapq.txt

    bedtools bamtobed -i ${genome_bam} \\
        | sort -k1,1 -k2,2n \\
        >  ${fasta_name}_genome.bed

    if [[ \$(wc -l  ${fasta_name}_genome.bed | awk '{print \$1}') -gt 0 ]]
    then
        cut -f4,6  ${fasta_name}_genome.bed \\
            | sort \\
            >  ${fasta_name}_genome_strand.txt
        cut -f4,5 ${fasta_name}_genome.bed \\
            | sort \\
            >  ${fasta_name}_genome_mapq.txt
        bedtools intersect -a  ${fasta_name}_genome.bed -b ${gene_bed} -wb \\
            > intersect_genome_genomes.bed

        if [[ \$(wc -l  intersect_genome_genomes.bed | awk '{print \$1}') -gt 0 ]]
        then
            cut -f1-4,10,6 intersect_genome_genomes.bed \\
                | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
                | cut -f4,6 \\
                | sort \\
                >  ${fasta_name}_genome_genes.txt
        fi
    fi

    genome_annotations.py \\
        --fasta_name ${fasta_name}
    """
}
