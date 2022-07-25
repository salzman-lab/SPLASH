
process GENOME_ANNOTATIONS {

    tag "${fasta_name}"
    label 'process_low'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6 pandas' : null)

    input:
    tuple path(fasta), path(end_to_end_genome_bam), path(end_to_end_transcriptome_bam), path(local_genome_bam), path(local_transcriptome_bam)
    path gene_bed

    output:
    path "*tsv"         , emit: annotations

    script:
    fasta_name          = fasta.baseName
    anchor_hits         = "genome_annotations_anchor.tsv"
    target_hits         = "genome_annotations_target.tsv"
    """
    ##
    ## really messy for now for diagnostic purposes
    ##


    ##
    ## local
    ##
    mkdir -p local
    samtools view ${local_transcriptome_bam} \\
        | cut -f1,3 \\
        > local/${fasta_name}_transcriptome_hit.txt

    samtools view ${local_transcriptome_bam} \\
        | cut -f1,5 \\
        > local/${fasta_name}_transcriptome_mapq.txt

    bedtools bamtobed -i ${local_genome_bam} \\
        | sort -k1,1 -k2,2n \\
        > local/${fasta_name}_genome.bed

    if [[ \$(wc -l local/${fasta_name}_genome.bed | awk '{print \$1}') -gt 0 ]]
    then
        cut -f4,6 local/${fasta_name}_genome.bed \\
            | sort \\
            > local/${fasta_name}_genome_strand.txt
        cut -f4,5 local/${fasta_name}_genome.bed \\
            | sort \\
            > local/${fasta_name}_genome_mapq.txt
        bedtools intersect -a local/${fasta_name}_genome.bed -b ${gene_bed} -wb \\
            | cut -f1-4,10,6 \\
            | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
            | cut -f4,6 \\
            | sort \\
            > local/${fasta_name}_genome_genes.txt
    fi

    ##
    ## end_to_end
    ##
    mkdir -p end_to_end
    samtools view ${end_to_end_transcriptome_bam} \\
        | cut -f1,3 \\
        > end_to_end/${fasta_name}_transcriptome_hit.txt

    samtools view ${end_to_end_transcriptome_bam} \\
        | cut -f1,5 \\
        >  end_to_end/${fasta_name}_transcriptome_mapq.txt

    bedtools bamtobed -i ${end_to_end_genome_bam} \\
        | sort -k1,1 -k2,2n \\
        >  end_to_end/${fasta_name}_genome.bed

    if [[ \$(wc -l  end_to_end/${fasta_name}_genome.bed | awk '{print \$1}') -gt 0 ]]
    then
        cut -f4,6  end_to_end/${fasta_name}_genome.bed \\
            | sort \\
            >  end_to_end/${fasta_name}_genome_strand.txt
        cut -f4,5 end_to_end/${fasta_name}_genome.bed \\
            | sort \\
            >  end_to_end/${fasta_name}_genome_mapq.txt
        bedtools intersect -a  end_to_end/${fasta_name}_genome.bed -b ${gene_bed} -wb \\
            | cut -f1-4,10,6 \\
            | bedtools groupby -i - -g 1,2,3,4,5 -c 6 -o collapse \\
            | cut -f4,6 \\
            | sort \\
            >  end_to_end/${fasta_name}_genome_genes.txt
    fi

    genome_annotations.py \\
        --fasta_name ${fasta_name}
    """
}
