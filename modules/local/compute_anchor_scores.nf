
process COMPUTE_ANCHOR_SCORES {

    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1" : null)


    input:
    path ch_targets_samplesheet
    val bound_distance
    val max_distance

    output:
    path outfile_counts_distances   , emit: anchor_target_counts
    path outfile_anchor_scores      , emit: anchor_scores
    path outfile_anchor_fasta       , emit: anchor_fasta
    path outfile_target_fasta       , emit: target_fasta

    script:
    outfile_counts_distances        = "anchor_targets_counts.tsv"
    outfile_anchor_scores           = "anchor_scores.tsv"
    outfile_anchor_fasta            = "anchors.fasta"
    outfile_target_fasta            = "targets.fasta"
    """
    compute_anchor_scores.py \\
        --samplesheet ${ch_targets_samplesheet} \\
        --bound_distance ${bound_distance} \\
        --max_distance ${max_distance} \\
        --outfile_counts_distances ${outfile_counts_distances} \\
        --outfile_anchor_scores ${outfile_anchor_scores} \\
        --outfile_anchor_fasta ${outfile_anchor_fasta} \\
        --outfile_target_fasta ${outfile_target_fasta}
    """
}
