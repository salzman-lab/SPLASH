
process GET_ABUNDANT_ANCHORS {

    label 'process_medium'

    input:
    path stratified_anchors
    val anchor_count_threshold
    val kmer_size

    output:
    path "*gz"                          , emit: abundant_stratified_anchors
    path outfile_abundant_control_seqs  , emit: abundant_control_seqs

    script:
    outfile_abundant_control_seqs        = "counts_abundant_${stratified_anchors}"
    """
    ## Create third column of the anchor sequence, so that the file = 54mer-sample counts, 54mer, anchor
    ## Sort for joining later
    paste -d' ' ${stratified_anchors} <(cut -f2 -d' ' ${stratified_anchors} | cut -c-${kmer_size}) \\
        | sort -k4 \\
        > with_anchors.txt

    ## Get file of 54mer-sample counts and anchors, and aggregate by summing over all anchors
    cut -f1,4 -d' ' with_anchors.txt \\
        | awk '{++freq[\$2]; sum[\$2]+=\$1} END{for (sample in sum) print sum[sample], sample}' \\
        > ${outfile_abundant_control_seqs}

    ## Filter sums by some threshold
    ## Sort of joining later
    awk -v anchor_count_threshold="${anchor_count_threshold}" '{if(\$1 > anchor_count_threshold) print \$2}' ${outfile_abundant_control_seqs} \\
        | sort \\
        > abundant_anchors.txt

    ## Join on abundant anchors
    ## Output same columns as input file
    join -1 1 -2 4 abundant_anchors.txt with_anchors.txt \\
        | cut -f2-4 -d' ' \\
        > abundant_${stratified_anchors}

    gzip -f abundant_${stratified_anchors}

    """
}
