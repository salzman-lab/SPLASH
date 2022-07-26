process UMITOOLS {

    tag "${R1.simpleName}"
    label "process_high"

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    // From nf-core
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38h4a8c8d9_0' :
        'quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0' }"

    input:
    tuple val(id), path(R1), path(R2)

    output:
    tuple val(id), path(extracted_R1), path(extracted_R2)   , emit: fastqs
    path "whitelist_${id}.txt"                              , emit: whitelist
    path "*log"                                             , emit: log

    script:
    extracted_R1 = "${id}_R1.extracted.fastq.gz"
    extracted_R2 = "${id}_R2.extracted.fastq.gz"
    """
    umi_tools whitelist \\
        --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNNNN \\
        --stdin ${R1} \\
        -L whitelist.${id}.log \\
        > whitelist_${id}.txt

    umi_tools extract \\
        --bc-pattern CCCCCCCCCCCCCCCCNNNNNNNNNNNN \\
        --stdin ${R1} \\
        --stdout ${extracted_R1} \\
        --read2-in ${R2} \\
        --read2-out ${extracted_R2} \\
        --whitelist=whitelist_${id}.txt
    """
}
