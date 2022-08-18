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
include { GET_FASTA                 } from '../modules/local/get_fasta'
include { GENOME_ALIGNMENT          } from '../modules/local/genome_alignment'
include { GENOME_ANNOTATIONS        } from '../modules/local/genome_annotations'
include { ELEMENT_ALIGNMENT         } from '../modules/local/element_alignment'
include { ELEMENT_ANNOTATIONS       } from '../modules/local/element_annotations'
include { SUMMARIZE_10X             } from '../modules/local/summarize_10X'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPROCESS_10X            } from '../subworkflows/local/preprocess_10X'
include { FETCH                     } from '../subworkflows/local/fetch'
include { PVALUES                   } from '../subworkflows/local/pvalues'


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


workflow NOMAD_10X {

    if (! params.abundant_stratified_anchors) {
        // Parse samplesheet
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

        PREPROCESS_10X(
            ch_paired_fastqs
        )

        // Define lookahead parameter
        if (params.lookahead) {

            lookahead = params.lookahead
        } else {

            // Use first fastq in sampelsheet to determine lookahead distance
            fastq = ch_paired_fastqs
                .first()
                .map{
                    it[2]
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
        // Subworkflow: Get abdunant anchors
        */
        FETCH(
            PREPROCESS_10X.out.fastqs,
            lookahead
        )

        abundant_control_seqs        = FETCH.out.abundant_control_seqs
        abundant_stratified_anchors  = FETCH.out.abundant_stratified_anchors

    } else {
        abundant_control_seqs        = Channel.empty()
        abundant_stratified_anchors  = Channel
            .fromPath("${params.abundant_stratified_anchors}/*")

    }

    /*
    // Subworkflow: Get control or significant anchors
    */
    PVALUES(
        abundant_control_seqs,
        abundant_stratified_anchors
    )

    // define
    anchors_pvals = PVALUES.out.anchors_pvals

    /*
    // Process: Make fasta from anchors
    */
    GET_FASTA(
        anchors_pvals
    )

    // Only proceed with anchor fasta file
    GET_FASTA.out.fasta
        .flatten()
        .filter{
            file ->
            file.name.contains('anchor.fasta')
        }
        .set { anchor_fasta }

    /*
    // Process: Align anchors to genome
    */
    GENOME_ALIGNMENT(
        anchor_fasta,
        params.genome_index,
        params.transcriptome_index
    )

    /*
    // Process: Run gene and exon annotations
    */
    GENOME_ANNOTATIONS(
        GENOME_ALIGNMENT.out.bam_tuple,
        params.gene_bed
    )

    // Parse samplesheet of bowtie2 indices
    element_annotations_samplesheet = file(
        params.element_annotations_samplesheet,
        checkIfExists: true
    )
    ch_indices = Channel.fromPath(element_annotations_samplesheet)
        .splitCsv(
            header: false
        )
        .map{ row ->
            row[0]
        }

    // Cartesian join of anchor fasta and all bowtie2 indices
    ch_anchor_indices = anchor_fasta
        .combine(ch_indices)

    /*
    // Process: Align anchors to each bowtie2 index
    */
    ELEMENT_ALIGNMENT(
        ch_anchor_indices
    )

    /*
    // Process: Merge scores with hits
    */
    ELEMENT_ANNOTATIONS(
        ELEMENT_ALIGNMENT.out.hits.collect()
    )

    /*
    // Process: Make summary file
    */
    SUMMARIZE_10X(
        anchors_pvals,
        ELEMENT_ANNOTATIONS.out.annotated_anchors,
        GENOME_ANNOTATIONS.out.annotated_anchors
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

