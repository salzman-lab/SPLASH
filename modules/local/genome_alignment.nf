
process GENOME_ALIGNMENT {

    tag "${index_name}"
    label 'process_low'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)

    input:
    path anchor_fasta
    path target_fasta
    val genome_index
    val transcriptome_index

    output:
    path anchor_genome_bam  , emit: anchor_genome_bam
    path target_genome_bam  , emit: target_genome_bam
    path anchor_trans_bam   , emit: anchor_trans_bam
    path target_trans_bam   , emit: target_trans_bam

    script:
    anchor_genome_bam       = "anchor_genome.bam"
    target_genome_bam       = "target_genome.bam"
    anchor_trans_bam        = "anchor_transcriptome.bam"
    target_trans_bam        = "target_transcriptome.bam"
    """
    bowtie2 -f -x ${genome_index} -U ${anchor_fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${anchor_genome_bam}

    bowtie2 -f -x ${genome_index} -U ${target_fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${target_genome_bam}

    bowtie2 -f -x ${transcriptome_index} -U ${anchor_fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${anchor_trans_bam}

    bowtie2 -f -x ${transcriptome_index} -U ${target_fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${target_trans_bam}
    """
}
