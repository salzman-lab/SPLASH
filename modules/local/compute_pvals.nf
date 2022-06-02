
process COMPUTE_PVALS {

    label 'process_high'

    input:
    path counts
    val kmer_size
    path samplesheet
    val K_num_hashes
    val L_num_random_Cj

    output:
    path outfile_scores     , emit: scores      , optional: true
    path "*extra_info*"     , emit: extra_info  , optional: true
    path "*pkl"             , emit: pkl         , optional: true


    script:
    file_id                 = counts.simpleName
    outfile_scores          = "scores_${file_id}.tsv"
    """
    compute_pvals.py \\
        --infile ${counts} \\
        --kmer_size ${kmer_size} \\
        --samplesheet ${samplesheet} \\
        --K ${K_num_hashes} \\
        --L ${L_num_random_Cj} \\
        --outfile ${outfile_scores}
    """
}
