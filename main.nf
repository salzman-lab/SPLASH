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
include { NOMAD_10X         } from './workflows/nomad_10X'
include { PLOT_HEATMAPS     } from './workflows/plot_heatmaps'

include { SAMPLESHEET_CHECK } from './modules/local/samplesheet_check'

//
// WORKFLOW: Run main kaitlinchaung/nomad analysis pipeline
//

workflow RUN_NOMAD {
    if (params.run_anchor_heatmaps){
        // Workflow to only plot heatmaps on previously-run data
        PLOT_HEATMAPS ()

    } else {

        if (params.is_10X) {
            // Workflow to run NOMAD on 10X samples
            NOMAD_10X ()

        } else {
            // Workflow to run NOMAD on bulk RNAseq and SS2 samples
            NOMAD ()

        }

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
