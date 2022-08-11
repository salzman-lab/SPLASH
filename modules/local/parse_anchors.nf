
process PARSE_ANCHORS {

    tag "${id}"
    label 'process_medium'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 conda-forge::biopython" : null)

    input:
    tuple val(id), path(fastq)
    path anchors_pvals_file
    val num_reads_second_pass
    val consensus_length
    val kmer_size
    val direction
    val lookahead

    output:
    path "*_target_counts.tsv"  , emit: targets
    path "*.fasta"              , emit: consensus_fasta
    path "*_counts.tab"         , emit: consensus_counts
    path "*_fractions.tab"      , emit: consensus_fractions
    path "*log"                 , emit: log

    script:
    out_consensus_fasta_file    = "${id}_${direction}.fasta"
    out_counts_file             = "${id}_${direction}_counts.tab"
    out_fractions_file          = "${id}_${direction}_fractions.tab"
    out_target_file             = "${id}_target_counts.tsv"
    """
    parse_anchors.py \\
        --num_lines ${num_reads_second_pass} \\
        --anchors_pvals_file ${anchors_pvals_file} \\
        --id ${id} \\
        --fastq_file ${fastq} \\
        --out_consensus_fasta_file ${out_consensus_fasta_file} \\
        --out_counts_file ${out_counts_file} \\
        --out_fractions_file ${out_fractions_file} \\
        --out_target_file ${out_target_file} \\
        --consensus_length ${consensus_length} \\
        --kmer_size ${kmer_size} \\
        --direction ${direction} \\
        --lookahead ${lookahead}
    """
}
