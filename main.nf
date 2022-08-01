#!/usr/bin/env nextflow
/*
========================================================================================
    kaitlinchaung/nomad
========================================================================================
    Github : https://github.com/kaitlinchaung/nomad

----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/


/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { NOMAD             } from './workflows/nomad'
include { ANCHOR_HEATMAPS   } from './modules/local/anchor_heatmaps'

//
// WORKFLOW: Run main kaitlinchaung/nomad analysis pipeline
//

workflow RUN_NOMAD {
    // If we are only plotting,
    if (params.run_anchor_heatmaps){

        // Declare file paths
        anchors_pvals                       = file("${params.results_dir}/anchors_pvals.tsv", checkIfExists: true)
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
            params.num_heatmap_anchors,
            file(params.input),
            additional_summary,
            genome_annotations_anchors
        )
    // Run pipeline
    } else {
        NOMAD ()

    }
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    RUN_NOMAD ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
