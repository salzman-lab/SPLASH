
process COMPUTE_ANCHOR_SCORES {

    publishDir "${outdir}", mode: 'copy'

    input:
    path ch_targets_samplesheet

    output:
    path outfile_counts_distances   , emit: anchor_target_counts
    path outfile_anchor_scores      , emit: anchor_scores
    path outfile_anchor_fasta       , emit: fasta

    script:
    outfile_counts_distances = "anchor_targets_counts.tsv"
    outfile_anchor_scores = "anchor_scores.tsv"
    outfile_anchor_fasta = "anchors.fasta"
    """
    compute_anchor_scores.py \\
        --samplesheet ${ch_targets_samplesheet} \\
        --outfile_counts_distances ${outfile_counts_distances} \\
        --outfile_anchor_scores ${outfile_anchor_scores} \\
        --outfile_anchor_fasta ${outfile_anchor_fasta}
    """
}
