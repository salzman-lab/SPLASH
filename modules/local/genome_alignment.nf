
process GENOME_ALIGNMENT {

    tag "${fasta_name}"
    label 'process_low'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)

    input:
    path fasta
    val genome_index
    val transcriptome_index

    output:
    tuple path(fasta), path(end_to_end_genome_bam), path(end_to_end_transcriptome_bam), path(local_genome_bam), path(local_transcriptome_bam), emit: bam_tuple

    script:
    fasta_name                      = fasta.baseName
    end_to_end_genome_bam           = "${fasta_name}_end_to_end_genome.bam"
    end_to_end_transcriptome_bam    = "${fasta_name}_end_to_end_transcriptome.bam"
    local_genome_bam                = "${fasta_name}_local_genome.bam"
    local_transcriptome_bam         = "${fasta_name}_local_transcriptome.bam"
    """
    bowtie2 -f -x ${genome_index} -U ${fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${end_to_end_genome_bam}

    bowtie2 -f -x ${genome_index} -U ${fasta} -k 1 --local --quiet \\
        | samtools view -bS - \\
        > ${local_genome_bam}

    bowtie2 -f -x ${transcriptome_index} -U ${fasta} -k 1 --quiet \\
        | samtools view -bS - \\
        > ${end_to_end_transcriptome_bam}

    bowtie2 -f -x ${transcriptome_index} -U ${fasta} -k 1 --local --quiet \\
        | samtools view -bS - \\
        > ${local_transcriptome_bam}
    """
}
