
process GENOME_ALIGNMENT {

    tag "${fasta_name}"
    label 'process_medium'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)

    input:
    path fasta
    val genome_index
    val transcriptome_index

    output:
    tuple path(fasta), path(genome_bam), path(transcriptome_bam), emit: bam_tuple
    path "*bam*"                    , emit: bam

    script:
    fasta_name                      = fasta.baseName
    genome_bam                      = "${fasta_name}_genome.bam"
    transcriptome_bam               = "${fasta_name}_transcriptome.bam"

    """
    bowtie2 -f -x ${genome_index} -U ${fasta} -k 1 -p ${task.cpus} --quiet \\
        | samtools view -bS - \\
        | samtools sort - \\
        > ${genome_bam}

    bowtie2 -f -x ${transcriptome_index} -U ${fasta} -k 1 -p ${task.cpus} --quiet \\
        | samtools view -bS - \\
        | samtools sort - \\
        > ${transcriptome_bam}

    for file in *bam
    do
        samtools index \${file}
    done
    """
}
