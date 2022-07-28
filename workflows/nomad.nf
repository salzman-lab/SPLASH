/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNomad.initialise(params, log)

// Check input path parameters to see if the files exist if they have been specified
checkPathParamList = [
    params.star_index,
    params.genome_index,
    params.gene_bed,
    params.gtf
]
for (param in checkPathParamList) {
    if (param) {
        file(param, checkIfExists: true)
        println(param)
    }
}

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
include { GET_READ_LENGTH           } from '../modules/local/get_read_length'
include { ADD_DUMMY_SCORE           } from '../modules/local/add_dummy_score'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FETCH_10X                 } from '../subworkflows/local/fetch_10X'
include { FETCH                     } from '../subworkflows/local/fetch'
include { ANALYZE                   } from '../subworkflows/local/analyze'
include { ANNOTATE                  } from '../subworkflows/local/annotate'


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

    // Validate samplesheet
    SAMPLESHEET_CHECK(
        file(params.input),
        params.is_10X
    )

    if (params.is_10X) {
        Channel
            .fromPath(params.input)
            .splitCsv(
                header: false
            )
            .map { row ->
                tuple(
                    row[0],
                    file(row[1]),
                    file(row[2])
                )
            }
            .set{ ch_paired_fastqs }

        FETCH_10X(
            ch_paired_fastqs
        )
        ch_fastqs = FETCH_10X.out.fastqs

    } else {
        // Parse samplesheet
        Channel
            .fromPath(params.input)
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
    }

    // Define lookahead parameter
    if (params.use_read_length) {
        /*
        // Get read length of dataset
        */
        if (params.is_10X){
            fastq = ch_paired_fastqs.first().map{
                it[2]
            }
        } else {
            fastq = ch_fastqs.first().map{
                it[1]
            }
        }
        GET_READ_LENGTH(
            fastq,
            SAMPLESHEET_CHECK.out.samplesheet
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

    if (! params.run_pvals_only) {
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
            ANALYZE.out.ch_consensus_fasta,
            ANALYZE.out.ch_anchor_target_fastas
        )
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
