include { GET_CONTROL_ANCHORS       } from '../../modules/local/get_control_anchors'
include { COMPUTE_PVALS             } from '../../modules/local/compute_pvals'
include { SIGNIFICANT_ANCHORS       } from '../../modules/local/significant_anchors'

workflow PVALUES {

    take:
    abundant_control_seqs
    abundant_stratified_anchors

    main:

    if (params.run_control) {
        GET_CONTROL_ANCHORS(
            abundant_control_seqs,
            params.num_control_anchors
        )

        anchors_pvals   = GET_CONTROL_ANCHORS.out.seqs
        anchors_Cjs     = Channel.empty()


    } else {
        /*
        // Process: Compute statistics
        */
        COMPUTE_PVALS(
            abundant_stratified_anchors,
            params.is_10X,
            params.kmer_size,
            file(params.cell_barcode_samplesheet),
            params.K_num_hashes,
            params.L_num_random_Cj,
            params.anchor_count_threshold,
            params.anchor_unique_targets_threshold,
            params.anchor_samples_threshold,
            params.anchor_sample_counts_threshold,
            params.anchor_batch_size
        )

        /*
        // Process: Output significant anchors
        */
        SIGNIFICANT_ANCHORS(
            COMPUTE_PVALS.out.scores.collect(),
            params.fdr_threshold,
            file(params.input)
        )

        // Only proceed with non empty files
        anchors_pvals   = SIGNIFICANT_ANCHORS.out.scores
            .filter{
                it.countLines() > 1
            }

        anchors_Cjs     = SIGNIFICANT_ANCHORS.out.cjs

    }

    emit:
    anchors_pvals               = anchors_pvals
    anchors_Cjs                 = anchors_Cjs

}
