/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNomad.initialise(params, log)


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { ANCHOR_HEATMAPS } from '../modules/local/anchor_heatmaps'


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/


workflow PLOT_HEATMAPS {

    // Declare file paths
    anchors_pvals                       = file("${params.results_dir}/anchors_pvals.tsv", checkIfExists: true)
    anchors_Cjs                         = file("${params.results_dir}/anchors_Cjs_random_opt.tsv", checkIfExists: true)
    genome_annotations_anchors          = "${params.results_dir}/genome_annotations/genome_annotations_anchor.tsv"
    additional_summary                  = "${params.results_dir}/additional_summary.tsv"
    abundant_stratified_anchors_path    = "${params.results_dir}/abundant_stratified_anchors/*txt.gz"
    consensus_fractions_path            = "${params.results_dir}/consensus_anchors/*fractions.tab"

    // Create channels with files
    Channel
        .fromPath(abundant_stratified_anchors_path)
        .set{ abundant_stratified_anchors}
    Channel
        .fromPath(consensus_fractions_path)
        .set{ consensus_fractions }

    if (params.use_heatmap_anchor_list) {
        heatmap_anchor_list = file(params.heatmap_anchor_list)
    } else {
        // Set a dummy value if we are not using heatmap anchor list
        heatmap_anchor_list = consensus_fractions.first().collectFile(name: 'dummy.txt')
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
        anchors_Cjs,
        params.num_heatmap_anchors,
        file(params.input),
        additional_summary,
        genome_annotations_anchors
    )

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
