process TRIMGALORE {
    tag "${fastq_id}"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::trim-galore=0.6.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0' :
        'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0' }"

    input:
    tuple val(fastq_id), path(fastq)

    output:
    tuple val(fastq_id), path("*trimmed.fq.gz") , emit: fastq
    tuple val(fastq_id), path("*report.txt")    , emit: log

    tuple val(fastq_id), path("*.html")         , emit: html optional true
    tuple val(fastq_id), path("*.zip")          , emit: zip optional true

    when:
    task.ext.when == null || task.ext.when

    script:

    // Clipping presets have to be evaluated in the context of SE/PE
    def c_r1   = params.clip_r1 > 0             ? "--clip_r1 ${params.clip_r1}"                         : ''
    def c_r2   = params.clip_r2 > 0             ? "--clip_r2 ${params.clip_r2}"                         : ''
    def tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    def tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''

    // Added soft-links to original fastqs for consistent naming in MultiQC
    """
    trim_galore \\
        ${fastq} \\
        --gzip
    """
}