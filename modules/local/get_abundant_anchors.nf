
process GET_ABUNDANT_ANCHORS {

    label 'process_medium'

    input:
    path counts
    val anchor_count_threshold

    output:
    path outfile_anchors            , emit: seqs
    path outfile_anchors_counts     , emit: anchor_counts

    script:
    outfile_anchors                 = "abundant_${counts}.gz"
    outfile_anchors_counts          = "counts_abundant_${counts}"
    """
    ## Create third column of the anchor sequence, so that the file = 54mer-sample counts, 54mer, anchor
    ## Sort for joining later
    paste -d' ' ${counts} <(cut -f2 -d' ' ${counts} | cut -c-27) \\
        | sort -k4 \\
        > with_anchors.txt

    ## Get file of 54mer-sample counts and anchors, and aggregate by summing over all anchors
    cut -f1,4 -d' ' with_anchors.txt \\
        | awk '{++freq[\$2]; sum[\$2]+=\$1} END{for (sample in sum) print sum[sample], sample}' \\
        > ${outfile_anchors_counts}

    ## Filter sums by some threshold
    ## Sort of joining later
        awk -v anchor_count_threshold="${anchor_count_threshold}" '{if(\$1 > anchor_count_threshold) print \$2}' ${outfile_anchors_counts} \\
        | sort \\
        > abundant_anchors.txt

    ## Join on abundant anchors
    ## Output same columns as input file
    join -1 1 -2 4 abundant_anchors.txt with_anchors.txt \\
        | cut -f2-4 -d' ' \\
        > ${outfile_anchors}

    gzip ${outfile_anchors}

    """
}
