include { PARSE_ANCHORS     } from '../../modules/local/parse_anchors'
include { ANCHORS_TARGETS   } from '../../modules/local/anchors_targets'


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

    /*
    // Process to get anchor scores and anchor-target counts
    */
    ANCHORS_TARGETS(
        PARSE_ANCHORS.out.targets.collect()
    )

    emit:
    anchors_targets         = ANCHORS_TARGETS.out.anchors_targets
    ch_consensus_fasta      = PARSE_ANCHORS.out.consensus_fasta.collect()
    ch_anchor_target_fastas = ANCHORS_TARGETS.out.fasta

}
