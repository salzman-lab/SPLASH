
process SUMMARIZE_10X {

    tag "${samplesheet_id}"
    label 'process_high_memory'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 numpy=1.22.3 bioconda::blast=2.12.2 bioconda::biopython=1.70" : null)

    input:
    tuple val(samplesheet_id), path(samplesheet), path(element_annotations), path(genome_annotations), path(anchors_pvals)

    output:
    path outfile        , emit: tsv

    script:
    outfile             = "summary_${samplesheet_id}.tsv"
    """
    summarize_10X.py \\
        --anchors_pvals ${anchors_pvals} \\
        --genome_annotations ${genome_annotations} \\
        --element_annotations ${element_annotations} \\
        --outfile ${outfile}
    """
}
