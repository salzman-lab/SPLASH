include { ANCHOR_HEATMAPS   } from '../../modules/local/anchor_heatmaps'


workflow PLOT {

    take:
    abundant_stratified_anchors
    consensus_fractions
    anchors_pvals
    additional_summary
    genome_annotations_anchors

    main:

    if (params.use_heatmap_anchor_list) {
        heatmap_anchor_list = file(params.heatmap_anchor_list)
    } else {
        // Set a dummy value if we are not using heatmap anchor list
        heatmap_anchor_list = consensus_fractions.collectFile(name: 'dummy.txt')
    }

    /*
    // Process to make heatmaps
    */
    ANCHOR_HEATMAPS(
        params.use_heatmap_anchor_list,
        heatmap_anchor_list,
        abundant_stratified_anchors.collect(),
        consensus_fractions.collect(),
        params.dataset,
        anchors_pvals,
        params.num_heatmap_anchors,
        file(params.input),
        additional_summary,
        genome_annotations_anchors
    )

}
