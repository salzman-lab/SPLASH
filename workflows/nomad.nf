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
include { GET_SOFTWARE_VERSIONS     } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { SAMPLESHEET_CHECK         } from '../modules/local/samplesheet_check'
include { GET_LOOKAHEAD             } from '../modules/local/get_lookahead'
include { ADD_DUMMY_SCORE           } from '../modules/local/add_dummy_score'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FETCH                     } from '../subworkflows/local/fetch'
include { PVALUES                   } from '../subworkflows/local/pvalues'
include { ANALYZE                   } from '../subworkflows/local/analyze'
include { ANNOTATE                  } from '../subworkflows/local/annotate'
include { PLOT                      } from '../subworkflows/local/plot'


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


workflow NOMAD {

    // Parse samplesheet
    Channel
        .fromPath(file(params.input))
        .splitCsv(
            header: false
        )
        .map { row ->
            tuple(
                file(row[0]).simpleName,
                file(row[0])
            )
        }
        .set{ ch_fastqs }


    // Define lookahead parameter
    if (params.lookahead) {
        lookahead = params.lookahead

    } else {
        // Use first fastq in sampelsheet to determine lookahead distance
        fastq = ch_fastqs
            .first()
            .map{
                it[1]
            }

        /*
        // Process: Get read lookahead distance of dataset
        */
        GET_LOOKAHEAD(
            fastq,
            file(params.input),
            params.kmer_size
        )

        lookahead = GET_LOOKAHEAD.out.lookahead.toInteger()
    }

    /*
    // Subworkflow: Get abundant anchors
    */
    FETCH(
        ch_fastqs,
        lookahead
    )

    /*
    // Subworkflow: Get control anchors
    */
    PVALUES(
        file(params.input),
        FETCH.out.abundant_control_seqs,
        FETCH.out.abundant_stratified_anchors
    )

    if (params.anchors_file || params.run_control) {

        if (params.run_control) {
            // If decoy, use control pvalues
            anchors = PVALUES.out.anchors_pvals

        } else {
            // If using anchors file, use file from input
            anchors = file(params.anchors_file)

        }

        /*
        // Process: Add a dummy scores columns for postprocessing
        */
        ADD_DUMMY_SCORE(
            anchors
        )

        anchors_pvals = ADD_DUMMY_SCORE.out.anchors_pvals

    } else {
        // Proceed with computed significant anchors
        anchors_pvals   = PVALUES.out.anchors_pvals
        anchors_Cjs     = PVALUES.out.anchors_Cjs

    }

    // Only proceed if we are doing more than pval caluclations
    if (! params.run_pvals_only) {
        /*
        // Subworkflow: Perform analysis
        */
        ANALYZE(
            ch_fastqs,
            anchors_pvals,
            lookahead
        )

        /*
        // Subworkflow: Perform annotations
        */
        ANNOTATE(
            anchors_pvals,
            ANALYZE.out.anchor_target_counts,
            ANALYZE.out.ch_consensus_fasta,
            ANALYZE.out.ch_anchor_target_fastas
        )

        if (params.anchors_file == null & params.run_control == false) {
            abundant_stratified_anchors = FETCH.out.abundant_stratified_anchors
            consensus_fractions         = ANALYZE.out.consensus_fractions
            additional_summary          = ANNOTATE.out.additional_summary
            genome_annotations_anchors  = ANNOTATE.out.genome_annotations_anchors

            // If annotations are run OR we only want to plot, run plot
            if (params.run_annotations){
                /*
                // Subworkflow: Perform plotting
                */
                PLOT(
                    abundant_stratified_anchors,
                    consensus_fractions,
                    anchors_pvals,
                    anchors_Cjs,
                    additional_summary,
                    genome_annotations_anchors
                )
            }
        }
    }
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
