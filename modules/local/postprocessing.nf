
process POSTPROCESSING {

    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 numpy=1.22.3 bioconda::blast=2.12.2 bioconda::biopython=1.70" : null)

    input:
    path scores
    path counts
    path annotated_anchors
    path annotated_targets
    val run_blast

    output:
    path ("*tsv")       , emit: tsv

    script:
    def run_blast = run_blast ? "--run_blast" : ""

    scores_hamming      = "anchor_norm_scores_hamming.tsv"
    counts_hamming      = "anchor_targets_counts_hamming.tsv"
    outfile_hamming     = "summary_annotations_hamming.tsv"

    scores_jaccard      = "anchor_norm_scores_jaccard.tsv"
    counts_jaccard      = "anchor_targets_counts_jaccard.tsv"
    outfile_jaccard     = "summary_annotations_jaccard.tsv"
    """
    postprocessing.py \\
        --anchor_norm_scores ${scores_hamming} \\
        --anchor_target_counts ${counts_hamming} \\
        --annotated_anchors ${annotated_anchors} \\
        --annotated_targets ${annotated_targets} \\
        --outfile ${outfile_hamming} \\
        ${run_blast}

    postprocessing.py \\
        --anchor_norm_scores ${scores_jaccard} \\
        --anchor_target_counts ${counts_jaccard} \\
        --annotated_anchors ${annotated_anchors} \\
        --annotated_targets ${annotated_targets} \\
        --outfile ${outfile_jaccard} \\
        ${run_blast}
    """
}
