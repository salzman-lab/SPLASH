
include { GET_READ_LENGTH           } from '../../modules/local/get_read_length'
include { SPLIT_FASTQS              } from '../../modules/local/split_fastqs'
include { GET_ANCHORS               } from '../../modules/local/get_anchors'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { COMPUTE_ANCHOR_SCORES     } from '../../modules/local/compute_anchor_scores'
include { BOWTIE2_ANNOTATION        } from '../../modules/local/bowtie2_annotation'
include { MERGE_ANCHOR_HITS         } from '../../modules/local/merge_anchor_hits'

workflow ANALYZE_FASTQS {
    take:
    samplesheet

    main:

    // parse samplesheet
    Channel
        .fromPath(samplesheet)
        .splitCsv(
            header: false
        )
        .map { row ->
            tuple(
                row[0],
                file(row[1])
            )
        }
        .set{ ch_fastqs }

    // definitions
    num_lines = params.chunk_size * params.n_iterations * 4
    num_chunk_lines = params.chunk_size * 4

    if (params.use_read_len) {
        /*
        // Get read length of dataset
        */
        GET_READ_LENGTH(
            ch_fastqs.first()
        )
        read_length = GET_READ_LENGTH.out.read_length.toInteger()
        looklength = ((read_length - 2 * params.kmer_size) / 2).toInteger()

    } else {
        looklength = params.looklength
    }

    if (params.anchors_file) {
        // use input anchors file
        ch_anchors = params.anchors_file

    } else {
        /*
        // Process to split each fastq into size read_chunk and gzip compress
        */
        SPLIT_FASTQS(
            ch_fastqs,
            num_lines,
            num_chunk_lines
        )

        ch_split_fastqs = SPLIT_FASTQS.out.fastq.collect()

        /*
        // Process to get list of candidate anchors via onthefly
        */

        GET_ANCHORS(
            ch_split_fastqs,
            params.n_iterations,
            params.chunk_size,
            params.kmer_size,
            file(samplesheet),
            params.target_counts_threshold,
            params.anchor_counts_threshold,
            params.anchor_freeze_threshold,
            params.anchor_score_threshold,
            params.anchor_mode,
            params.c_type,
            params.window_slide,
            looklength,
            params.num_keep_anchors,
            params.use_std,
            params.compute_target_distance
        )

        ch_anchors = GET_ANCHORS.out.anchors
    }

    /*
    // Process to get consensus sequences and target counts for annchors
    */
    PARSE_ANCHORS(
        ch_fastqs,
        ch_anchors,
        num_lines,
        params.consensus_length,
        params.kmer_size,
        params.direction,
        looklength
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
    COMPUTE_ANCHOR_SCORES(
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
    // Process to align anchors to each bowtie2 index
    */
    BOWTIE2_ANNOTATION(
        COMPUTE_ANCHOR_SCORES.out.fasta.first(),
        ch_indices
    )

    // create samplesheet of the anchor hits files
    BOWTIE2_ANNOTATION.out.anchor_hits
        .collectFile() { file ->
            def X=file; X.toString() + '\n'
        }
        .set{ anchor_hits_samplesheet }

    /*
    // Process to merge anchor scores with anchor hits
    */
    MERGE_ANCHOR_HITS(
        anchor_hits_samplesheet,
        COMPUTE_ANCHOR_SCORES.out.anchor_scores
    )

}


