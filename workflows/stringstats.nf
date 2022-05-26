/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowStringstats.initialise(params, log)


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
include { GET_READ_LENGTH           } from '../modules/local/get_read_length'
include { ADD_DUMMY_SCORE           } from '../modules/local/add_dummy_score'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FETCH             } from '../subworkflows/local/fetch'
include { ANALYZE           } from '../subworkflows/local/analyze'
include { ANNOTATE          } from '../subworkflows/local/annotate'



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


workflow STRINGSTATS {

    // Parse samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(
            header: false
        )
        .map { row ->
            tuple(
                file(row[0]).simpleName,
                file(row[0]),
                row[1]
            )
        }
        .set{ ch_fastqs }

    // Define lookahead parameter
    if (params.use_read_length) {
        /*
        // Get read length of dataset
        */
        GET_READ_LENGTH(
            ch_fastqs.first()
        )
        read_length = GET_READ_LENGTH.out.read_length.toInteger()
        lookahead = ((read_length - 2 * params.kmer_size) / 2).toInteger()

    } else {
        lookahead = params.lookahead
    }

    if (params.anchors_file || params.run_decoy) {

        if (params.run_decoy) {
            /*
            // Get abundant anchors
            */
            FETCH(
                ch_fastqs,
                lookahead
            )

            anchors = FETCH.out.anchors_scores

        } else {
            anchors = file(params.anchors_file)
        }

        /*
        // Add a dummy scores columns for postprocessing
        */
        ADD_DUMMY_SCORE(
            anchors
        )

        anchors_scores = ADD_DUMMY_SCORE.out.anchors_scores

    } else {
        /*
        // Get significant anchors
        */
        FETCH(
            ch_fastqs,
            lookahead
        )

        anchors_scores = FETCH.out.anchors_scores

    }

    /*
    // Perform analysis
    */
    ANALYZE(
        ch_fastqs,
        anchors_scores,
        lookahead
    )

    /*
    // Perform annotations
    */
    ANNOTATE(
        anchors_scores,
        ANALYZE.out.anchor_target_counts,
        ANALYZE.out.ch_consensus_fasta
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
