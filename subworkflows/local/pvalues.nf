include { GET_CONTROL_ANCHORS       } from '../../modules/local/get_control_anchors'
include { MAKE_SAMPLESHEETS         } from '../../modules/local/make_samplesheets'
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
        // Process: Make pairwise samplesheets for pvalue computation
        */
        if (params.is_10X) {
            base_samplesheet = file(params.cell_barcode_samplesheet)

        } else {
            base_samplesheet = file(params.input)

        }

        MAKE_SAMPLESHEETS(
            base_samplesheet,
            params.is_10X,
            params.run_unsupervised_pvals
        )

        // Map a samplesheet id to each samplesheet
        ch_pairwise_samplesheets = MAKE_SAMPLESHEETS.out.samplesheets
            .flatten()
            .map { it ->
                if (params.run_unsupervised_pvals) {
                    tuple(
                        "unsupervised_pvalues",
                        it
                    )
                } else {
                    tuple(
                        it.simpleName.replaceFirst(/pairwise_samplesheet_/, ''),
                        it
                    )
                }
            }
            .view()

        // Carteisian product of all pairwise samplesheets and abundant stratified anchors files
        ch_samplesheets_anchors = ch_pairwise_samplesheets
            .combine(abundant_stratified_anchors)

        /*
        // Process: Compute statistics
        */
        COMPUTE_PVALS(
            ch_samplesheets_anchors,
            params.run_unsupervised_pvals,
            params.kmer_size,
            params.K_num_hashes,
            params.L_num_random_Cj,
            params.anchor_count_threshold,
            params.anchor_unique_targets_threshold,
            params.anchor_samples_threshold,
            params.anchor_sample_counts_threshold,
            params.anchor_batch_size
        )

        // For each samplesheet id, map to samplesheet and all associated abudant stratified anchor files
        // We only need to map to the first samplesheet of the groupTuple because they are all the same file but just in diff workDirs
        ch_samplesheets_scores = COMPUTE_PVALS.out.scores
            .groupTuple()
            .map{ it ->
                tuple(
                    it[0],
                    it[1][0],
                    it[2]
                )
            }

        /*
        // Process: Output significant anchors
        */
        SIGNIFICANT_ANCHORS(
            ch_samplesheets_scores,
            params.fdr_threshold
        )

        // Only proceed with non empty pvals files
        anchors_Cjs   = SIGNIFICANT_ANCHORS.out.cjs
        anchors_pvals = SIGNIFICANT_ANCHORS.out.scores
            .filter{
                it[2].countLines() > 1
            }

    }

    emit:
    anchors_pvals               = anchors_pvals
    anchors_Cjs                 = anchors_Cjs

}
