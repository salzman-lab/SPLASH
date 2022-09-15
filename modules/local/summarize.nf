
process SUMMARIZE {

    tag "${samplesheet_id}"
    label 'process_high_memory'
    publishDir(
        path: {"${params.outdir}/${samplesheet_id}"},
        mode: "copy",
        pattern: "*.tsv")
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 numpy=1.22.3 bioconda::blast=2.12.2 bioconda::biopython=1.70" : null)

    input:
    tuple val(samplesheet_id), path(anchors_pvals), path(element_annotations), path(genome_annotations)

    output:
    path outfile        , emit: tsv

    script:
    outfile             = "summary.tsv"
    """
    summarize.py \\
        --anchors_pvals ${anchors_pvals} \\
        --genome_annotations ${genome_annotations} \\
        --element_annotations ${element_annotations} \\
        --outfile ${outfile}
    """
}
