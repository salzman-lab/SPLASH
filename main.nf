#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/stringstats
========================================================================================
    Github : https://github.com/nf-core/stringstats
    Website: https://nf-co.re/stringstats
    Slack  : https://nfcore.slack.com/channels/stringstats
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

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

include { STRINGSTATS } from './workflows/stringstats'

//
// WORKFLOW: Run main nf-core/stringstats analysis pipeline
//
workflow NFCORE_STRINGSTATS {
    STRINGSTATS ()
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
    NFCORE_STRINGSTATS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
