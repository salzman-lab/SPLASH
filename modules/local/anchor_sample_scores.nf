
process ANCHOR_SAMPLE_SCORES {

    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1" : null)

    input:
    path targets_samplesheet
    val bound_distance
    val max_distance
    val kmer_size
    val distance_type

    output:
    path outfile_counts_distances       , emit: anchor_target_counts
    path outfile_anchor_sample_scores   , emit: anchor_sample_scores

    script:
    outfile_counts_distances               = "anchor_targets_counts.tsv"
    outfile_anchor_sample_scores           = "anchor_sample_scores.tsv"
    """
    anchor_sample_scores.py \\
        --targets_samplesheet ${targets_samplesheet} \\
        --bound_distance ${bound_distance} \\
        --max_distance ${max_distance} \\
        --kmer_size ${kmer_size} \\
        --distance_type ${distance_type} \\
        --outfile_counts_distances ${outfile_counts_distances} \\
        --outfile_anchor_sample_scores ${outfile_anchor_sample_scores}
    """
}
