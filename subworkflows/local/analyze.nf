include { PARSE_ANCHORS         } from '../../modules/local/parse_anchors'
include { ANCHOR_TARGET_COUNTS  } from '../../modules/local/anchor_target_counts'


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

    // Create samplesheet of target counts files
    PARSE_ANCHORS.out.targets
        .collectFile() { file ->
            def X=file; X.toString() + '\n'
        }
        .set{ targets_samplesheet }

    /*
    // Process to get anchor scores and anchor-target counts
    */
    ANCHOR_TARGET_COUNTS(
        targets_samplesheet
    )


    emit:
    anchor_target_counts    = ANCHOR_TARGET_COUNTS.out.anchor_target_counts.first()
    ch_consensus_fasta      = PARSE_ANCHORS.out.consensus_fasta.collect()

}
