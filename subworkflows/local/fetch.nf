include { FETCH_ANCHORS             } from '../../modules/local/fetch_anchors'
include { COUNT_ANCHORS             } from '../../modules/local/count_anchors'
include { STRATIFY_ANCHORS          } from '../../modules/local/stratify_anchors'
include { GET_ABUNDANT_ANCHORS      } from '../../modules/local/get_abundant_anchors'

workflow FETCH {

    take:
    ch_fastqs
    lookahead

    main:

    /*
    // Process: Get all candidate anchors and targets
    */
    FETCH_ANCHORS(
        ch_fastqs,
        params.is_10X,
        params.num_reads_first_pass,
        params.kmer_size,
        lookahead,
        params.anchor_mode,
        params.window_slide
    )

    /*
    // Process: Count anchors and targets
    */
    COUNT_ANCHORS(
        FETCH_ANCHORS.out.seqs
    )

    /*
    // Process: Stratify counts files by 3mers
    */
    STRATIFY_ANCHORS(
        COUNT_ANCHORS.out.seqs.collect(),
        params.stratify_level
    )

    // Do not proceed with TTT file if in RNAseq mode
    if (params.is_RNAseq) {
        ch_stratified_anchors = STRATIFY_ANCHORS.out.seqs
            .flatten()
            .filter {
                it -> !file(it).name.contains("TTT")
            }

    } else {
        ch_stratified_anchors = STRATIFY_ANCHORS.out.seqs.flatten()

    }

    /*
    // Process: Filter kmer counts for abundant anchors
    */
    GET_ABUNDANT_ANCHORS(
        ch_stratified_anchors,
        params.anchor_count_threshold,
        params.kmer_size
    )

    emit:
    abundant_control_seqs           = GET_ABUNDANT_ANCHORS.out.abundant_control_seqs
    abundant_stratified_anchors     = GET_ABUNDANT_ANCHORS.out.abundant_stratified_anchors.flatten()

}
