
process GENOME_ALIGNMENT {

    tag "${index_name}"
    label 'process_low'
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)

    input:
    path anchor_fasta
    path target_fasta
    val index

    output:
    path "anchor_*bam"  , emit: anchor_bam
    path "target_*bam"  , emit: target_bam

    script:
    index_name          = file(index).baseName
    anchor_bam          = "anchor_genome.bam"
    target_bam          = "target_genome.bam"
    """

    bowtie2 -f -x ${index} -U ${anchor_fasta} -k 1 --quiet \\
        | samtools view -bS -
        >> ${anchor_bam}

    bowtie2 -f -x ${index} -U ${target_fasta} -k 1 --quiet \\
        | samtools view -bS -
        >> ${target_bam}
    """
}
