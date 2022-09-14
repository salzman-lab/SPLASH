
process GENOME_ALIGNMENT {
    tag "${samplesheet_id}"
    label 'process_medium'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)

    input:
    tuple val(samplesheet_id), path(samplesheet), path(fasta)
    val genome_index
    val transcriptome_index

    output:
    tuple val(samplesheet_id), path(samplesheet), path(fasta), path(genome_bam),
        path(transcriptome_bam)     , emit: bams
    path "*bam*"                    , emit: bam

    script:
    genome_bam                      = "${samplesheet_id}_genome.bam"
    transcriptome_bam               = "${samplesheet_id}_transcriptome.bam"
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
