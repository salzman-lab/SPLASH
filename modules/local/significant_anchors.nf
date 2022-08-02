
process SIGNIFICANT_ANCHORS {

    label 'process_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 anaconda::statsmodels" : null)

    input:
    path anchors
    val fdr_threshold
    path samplesheet

    output:
    path outfile_scores         , emit: scores
    path "all*tsv"              , emit: all_scores  , optional: true
    path "cjs*"                 , emit: cjs         , optional: true

    script:
    outfile_scores              = "anchors_pvals.tsv"
    outfile_all_anchors_pvals   = "all_anchors_pvals.tsv"
    outfile_Cjs                 = "cjs_rand_opt.tsv"
    """
    significant_anchors.py \\
        --fdr_threshold ${fdr_threshold} \\
        --samplesheet ${samplesheet} \\
        --outfile_scores ${outfile_scores} \\
        --outfile_all_anchors_pvals ${outfile_all_anchors_pvals} \\
        --outfile_Cjs ${outfile_Cjs}
    """
}
