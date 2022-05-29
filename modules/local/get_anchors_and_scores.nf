
process GET_ANCHORS_AND_SCORES {

    label 'process_medium'

    input:
    path counts
    val kmer_size

    output:
    path outfile_scores     , emit: scores
    path "*extra_info*"     , emit: extra_info
    path "*pkl"     , emit: pkl


    script:
    file_id                 = counts.simpleName
    outfile_scores          = "scores_${file_id}.tsv"
    """
    randHash_parallel.py \\
        --infile ${counts} \\
        --kmer_size ${kmer_size} \\
        --outfile ${outfile_scores}
    """
}
