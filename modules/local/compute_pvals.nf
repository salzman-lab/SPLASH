
process COMPUTE_PVALS {

    tag "${file_id}"
    label 'process_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.3 numpy conda-forge::mmh3=3.0.0 anaconda::scipy=1.7.3 anaconda::nltk=3.7" : null)

    input:
    tuple val(samplesheet_id), path(samplesheet), path(counts)
    val is_10X
    val kmer_size
    val K_num_hashes
    val L_num_random_Cj
    val anchor_count_threshold
    val anchor_unique_targets_threshold
    val anchor_samples_threshold
    val anchor_sample_counts_threshold
    val anchor_batch_size

    output:
    tuple val(samplesheet_id), path(samplesheet), path(outfile_scores), emit: scores, optional: true
    path "*extra_info*"     , emit: extra_info  , optional: true
    path "*pkl"             , emit: pkl         , optional: true

    script:
    def is_10X              = (is_10X == true)  ? "--is_10X" : ""
    file_id                 = counts.simpleName
    outfile_scores          = "scores_${file_id}.tsv"
    """
    compute_pvals.py \\
        --infile ${counts} \\
        --kmer_size ${kmer_size} \\
        --samplesheet ${samplesheet} \\
        --K ${K_num_hashes} \\
        --L ${L_num_random_Cj} \\
        --anchor_count_threshold ${anchor_count_threshold} \\
        --anchor_unique_targets_threshold ${anchor_unique_targets_threshold} \\
        --anchor_samples_threshold ${anchor_samples_threshold} \\
        --anchor_sample_counts_threshold ${anchor_sample_counts_threshold} \\
        --outfile ${outfile_scores} \\
        --anchor_batch_size ${anchor_batch_size} \\
        ${is_10X}
    """
}
