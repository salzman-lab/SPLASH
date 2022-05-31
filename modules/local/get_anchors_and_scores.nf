
process GET_ANCHORS_AND_SCORES {

    label 'process_high'

    input:
    path counts
    val kmer_size
    path samplesheet

    output:
    path outfile_scores     , emit: scores      , optional: true
    path "*extra_info*"     , emit: extra_info  , optional: true
    path "*pkl"             , emit: pkl         , optional: true


    script:
    file_id                 = counts.simpleName
    outfile_scores          = "scores_${file_id}.tsv"
    """
    randHash_parallel.py \\
        --infile ${counts} \\
        --kmer_size ${kmer_size} \\
        --samplesheetIDs ${samplesheet} \\
        --outfile ${outfile_scores}
    """
}
