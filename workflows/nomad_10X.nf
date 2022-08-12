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
include { SAMPLESHEET_CHECK         } from '../modules/local/samplesheet_check'
include { GET_READ_LENGTH           } from '../modules/local/get_read_length'
include { GET_FASTA                 } from '../modules/local/get_fasta'
include { GENOME_ALIGNMENT          } from '../modules/local/genome_alignment'
include { GENOME_ANNOTATIONS        } from '../modules/local/genome_annotations'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPROCESS_10X            } from '../subworkflows/local/preprocess_10X'
include { FETCH                     } from '../subworkflows/local/fetch'
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

    profile = "Project : $workflow.profile"

    if (profile.contains('test')) {
        // Skip the samplesheet check if this is test run
        ch_samplesheet = file(params.input)

    } else {
        /*
        // Process to validate sampelsheet
        */
        SAMPLESHEET_CHECK(
            file(params.input),
            params.is_10X
        )

        ch_samplesheet = SAMPLESHEET_CHECK.out.samplesheet

    }

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

    /*
    // Subworkflow to preprocess 10X reads
    */
    PREPROCESS_10X(
        file(params.input),
        ch_paired_fastqs
    )
    ch_fastqs = PREPROCESS_10X.out.fastqs


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
            ch_samplesheet
        )
        read_length = GET_READ_LENGTH.out.read_length.toInteger()
        lookahead = ((read_length - 2 * params.kmer_size) / 2).toInteger()
    } else {
        lookahead = params.lookahead
    }

    /*
    // Subworkflow: Get significant anchors
    */
    FETCH(
        ch_fastqs,
        lookahead
    )

    /*
    // Process: Make fasta from anchors
    */
    GET_FASTA(
        FETCH.out.anchors_pvals
    )

    anchor_fasta = GET_FASTA.out.annotations
        .filter{
            file ->
            file.name.contains('anchor.fasta')
        }

    /*
    // Process to align targets and anchors to genome
    */
    GENOME_ALIGNMENT(
        anchor_fasta,
        params.genome_index,
        params.transcriptome_index
    )

    /*
    // Process to run gene and exon annotations
    */
    GENOME_ANNOTATIONS(
        GENOME_ALIGNMENT.out.bam_tuple,
        params.gene_bed
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
