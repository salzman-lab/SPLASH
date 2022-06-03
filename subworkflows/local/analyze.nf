include { PARSE_ANCHORS         } from '../../modules/local/parse_anchors'
include { ANCHOR_TARGET_COUNTS  } from '../../modules/local/anchor_target_counts'
include { ANCHORS_TARGETS       } from '../../modules/local/anchors_targets'


workflow ANALYZE {

    take:
    ch_fastqs
    anchors_scores
    lookahead

    main:

    /*
    // Process to get consensus sequences and target counts for annchors
    */
    PARSE_ANCHORS(
        ch_fastqs,
        anchors_scores,
        params.num_reads_second_pass,
        params.consensus_length,
        params.kmer_size,
        params.direction,
        lookahead
    )

    if (params.run_anchor_target_counts) {
        /*
        // Process to get anchor scores and anchor-target counts
        */
        ANCHOR_TARGET_COUNTS(
            PARSE_ANCHORS.out.targets.collect()
        )

        anchor_target_counts    = ANCHOR_TARGET_COUNTS.out.anchor_target_counts
        ch_anchor_target_fastas = ANCHOR_TARGET_COUNTS.out.fasta

    } else {
        /*
        // Process to get anchors and their targets
        */
        ANCHORS_TARGETS(
            PARSE_ANCHORS.out.targets.collect()
        )

        anchor_target_counts    = ANCHORS_TARGETS.out.anchors_targets
        ch_anchor_target_fastas = ANCHORS_TARGETS.out.fasta
    }

    emit:
    anchor_target_counts        = anchor_target_counts
    ch_consensus_fasta          = PARSE_ANCHORS.out.consensus_fasta.collect()
    ch_anchor_target_fastas     = ch_anchor_target_fastas

}
