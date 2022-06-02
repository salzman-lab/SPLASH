
process ANCHOR_TARGET_COUNTS {

    label 'process_memory'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1" : null)

    input:
    path target_counts

    output:
    path "anchor_targets_counts*.tsv"   , emit: anchor_target_counts
    path "*.fasta"                      , emit: fasta

    script:
    outfile_anchor_target_counts        = "anchor_targets_counts.tsv"
    outfile_anchor_fasta                = "anchor.fasta"
    outfile_target_fasta                = "target.fasta"
    """
    anchor_target_counts.py \\
        --outfile_anchor_target_counts ${outfile_anchor_target_counts} \\
        --outfile_anchor_fasta ${outfile_anchor_fasta} \\
        --outfile_target_fasta ${outfile_target_fasta}
    """
}
