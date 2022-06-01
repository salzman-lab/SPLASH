
process SUMMARIZE {

    label 'process_extra_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 numpy=1.22.3 bioconda::blast=2.12.2 bioconda::biopython=1.70" : null)


    input:
    path anchor_scores
    path anchors_targets
    path annotated_anchors
    path annotated_targets
    val run_blast

    output:
    path outfile                , emit: tsv

    script:
    def run_blast = run_blast   ? "--run_blast" : ""

    outfile                     = "summary.tsv"

    """
    summarize.py \\
        --anchor_scores ${anchor_scores} \\
        --anchors_targets ${anchors_targets} \\
        --annotated_anchors ${annotated_anchors} \\
        --annotated_targets ${annotated_targets} \\
        --outfile ${outfile} \\
        ${run_blast}
    """
}
