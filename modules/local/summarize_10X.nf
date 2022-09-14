
process SUMMARIZE_10X {

    tag "${samplesheet_id}"
    label 'process_high_memory'
    publishDir(
        path: {"${params.outdir}/${samplesheet_id}"},
        mode: "copy",
        pattern: "*.tsv")
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 numpy=1.22.3 bioconda::blast=2.12.2 bioconda::biopython=1.70" : null)

    input:
    tuple val(samplesheet_id), path(genome_annotations), path(samplesheet), path(element_annotations), path(anchors_pvals)

    output:
    path outfile        , emit: tsv

    script:
    outfile             = "summary.tsv"
    """
    summarize_10X.py \\
        --outfile ${outfile}
    """
}
