include { ANCHOR_HEATMAPS   } from '../../modules/local/anchor_heatmaps'


workflow PLOT {

    take:
    abundant_stratified_anchors
    consensus_fractions
    anchors_pvals
    anchors_Cjs
    additional_summary
    genome_annotations_anchors

    main:

    // Create a dummly anchor list file, to re-use ANCHOR_HEATMAPS process
    heatmap_anchor_list = anchors_pvals.collectFile(name: 'dummy.txt')

    /*
    // Process to make heatmaps
    */
    ANCHOR_HEATMAPS(
        params.use_heatmap_anchor_list,
        heatmap_anchor_list,
        abundant_stratified_anchors.collect(),
        consensus_fractions.collect(),
        params.dataset,
        params.kmer_size,
        anchors_pvals,
        anchors_Cjs,
        params.num_heatmap_anchors,
        file(params.input),
        additional_summary,
        genome_annotations_anchors,
        params.results_dir
    )

}
