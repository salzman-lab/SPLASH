include { GET_READ_LENGTH           } from '../../modules/local/get_read_length'
include { GET_UNMAPPED              } from '../../modules/local/get_unmapped'
include { FETCH_ANCHORS             } from '../../modules/local/fetch_anchors'
include { COUNT_ANCHORS             } from '../../modules/local/count_anchors'
include { STRATIFY_ANCHORS          } from '../../modules/local/stratify_anchors'
include { GET_ANCHORS_AND_SCORES    } from '../../modules/local/get_anchors_and_scores'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { MERGE_TARGET_COUNTS       } from '../../modules/local/merge_target_counts'
include { GET_FASTA                 } from '../../modules/local/get_fasta'
include { BOWTIE2_ANNOTATION        } from '../../modules/local/bowtie2_annotation'
include { MERGE_ANNOTATIONS         } from '../../modules/local/merge_annotations'
include { POSTPROCESSING            } from '../../modules/local/postprocessing'

workflow SKIP_GET_ANCHORS {

    main:

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

    // /*
    // // Trim fastqs
    // */
    // TRIMGALORE(
    //     ch_fastqs
    // )

    // /*
    // // Check if we are only using unmapped reads
    // */
    // if (params.unmapped) {
    //     /*
    //     // Get unmapped reads
    //     */
    //     GET_UNMAPPED(
    //         TRIMGALORE.out.fastq,
    //         params.index_bowtie
    //     )

    //     ch_fastqs = GET_UNMAPPED.out.fastq

    // } else {
    //     ch_fastqs = TRIMGALORE.out.fastq

    // }


    /*
    // Process to get consensus sequences and target counts for annchors
    */
    PARSE_ANCHORS(
        ch_fastqs,
        file(params.anchors_file),
        params.num_parse_anchors_reads,
        params.consensus_length,
        params.kmer_size,
        params.direction,
        lookahead
    )

    // Create samplesheet of target counts files
    PARSE_ANCHORS.out.targets
        .collectFile() { file ->
            def X=file; X.toString() + '\n'
        }
        .set{ targets_samplesheet }

    /*
    // Process to get anchor scores and anchor-target counts
    */
    MERGE_TARGET_COUNTS(
        targets_samplesheet
    )

    // create samplesheet of bowtie2 indices
    ch_indices = Channel.fromPath(params.bowtie2_samplesheet)
        .splitCsv(
            header: false
        )
        .map{ row ->
            row[0]
        }

    /*
    // Process to get anchor and target fastas
    */
    GET_FASTA(
        MERGE_TARGET_COUNTS.out.anchor_target_counts.first()
    )

    /*
    // Process to align anchors to each bowtie2 index
    */
    BOWTIE2_ANNOTATION(
        GET_FASTA.out.anchor_fasta,
        GET_FASTA.out.target_fasta,
        ch_indices
    )

    // create samplesheet of the anchor hits files
    BOWTIE2_ANNOTATION.out.anchor_hits
        .collectFile(name: "anchor_samplesheet.txt") { file ->
            def X=file; X.toString() + '\n'
        }
        .set{ anchor_hits_samplesheet }

    // create samplesheet of the target hits files
    BOWTIE2_ANNOTATION.out.target_hits
        .collectFile(name: "target_samplesheet.txt") { file ->
            def X=file; X.toString() + '\n'
        }
        .set{ target_hits_samplesheet }

    /*
    // Process to merge scores with hits
    */
    MERGE_ANNOTATIONS(
        anchor_hits_samplesheet,
        target_hits_samplesheet
    )


}
