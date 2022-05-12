
process STAR_ALIGN {

    label 'process_high_memory'
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.5a--0' :
        'quay.io/biocontainers/star:2.7.5a--0' }"

    input:
    path fasta
    path star_index
    path gtf

    output:
    tuple val(fasta_name), path("*SJ.out.tab")                  , emit: junctions
    path "${fasta_name}__STARgenome/sjdbList.fromGTF.out.tab"   , emit: gtf_junctions

    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    fasta_name  = fasta.baseName
    """
    STAR \\
        --runThreadN ${cores} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fasta} \\
        --twopassMode Basic \\
        --alignIntronMax 1000000 \\
        --outFileNamePrefix ${fasta_name}_ \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes All \\
        --chimOutType WithinBAM SoftClip Junctions \\
        --chimJunctionOverhangMin 10 \\
        --chimSegmentReadGapMax 0 \\
        --chimOutJunctionFormat 1 \\
        --chimSegmentMin 12 \\
        --chimScoreJunctionNonGTAG -4 \\
        --chimNonchimScoreDropMin 10 \\
        --quantMode GeneCounts \\
        --sjdbGTFfile ${gtf} \\
        --outReadsUnmapped Fastx

    rm -rf ${fasta_name}__STARgenome/*txt
    rm -rf ${fasta_name}__STARgenome/exon*
    rm -rf ${fasta_name}__STARgenome/gene*
    rm -rf ${fasta_name}__STARgenome/transcript*
    rm -rf ${fasta_name}__STARgenome/sjdbList.out.tab
    rm -rf ${fasta_name}__STARpass1
    """
}
