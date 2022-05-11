
process GENOME_ALIGNMENT {

    tag "${index_name}"
    label 'process_low'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)

    input:
    path fasta
    val genome_index
    val transcriptome_index

    output:
    tuple path(fasta), path(genome_bam), path(transcriptome_bam), emit: bam_tuple

    script:
    fasta_name      = fasta.baseName
    genome_bam      = "${fasta_name}_genome.bam"
    transcriptome_bam       = "${fasta_name}_transcriptome.bam"
    """
    bowtie2 -f -x ${genome_index} -U ${fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${genome_bam}

    bowtie2 -f -x ${genome_index} -U ${fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${transcriptome_bam}
    """
}
