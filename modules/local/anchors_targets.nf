
process ANCHORS_TARGETS {

    label 'process_high_memory'

    input:
    path target_counts

    output:
    path "anchors_targets.tsv"  , emit: anchors_targets
    path "*.fasta"              , emit: fasta

    script:
    outfile_anchors_targets     = "anchors_targets.tsv"
    outfile_anchor_fasta        = "anchor.fasta"
    outfile_target_fasta        = "target.fasta"
    """
    rm -rf all_anchor_targets.tsv

    for file in *_target_counts.tsv
    do
        tail -n +2 \${file} \\
            | cut -f1,2 \\
            >> all_anchor_targets.tsv
    done

    sort all_anchor_targets.tsv | uniq > ${outfile_anchors_targets}

    cut -f1 ${outfile_anchors_targets} \
        | sort \
        | uniq \
        | awk '{print ">"\$0"\\n"\$0}' \
        > ${outfile_anchor_fasta}

    cut -f2 ${outfile_anchors_targets} \
        | sort \
        | uniq \
        | awk '{print ">"\$0"\\n"\$0}' \
        > ${outfile_target_fasta}
    """
}
