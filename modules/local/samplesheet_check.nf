
process SAMPLESHEET_CHECK {
    tag "${samplesheet}"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path samplesheet
    val is_10X

    output:
    path samplesheet    , emit: samplesheet

    script: // This script is bundled with the pipeline, in kaitlinchaung/nomad/bin/
    def is_10X          = (is_10X == true) ? "--is_10X" : ""
    """
    check_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        ${is_10X}
    """
}
