include { GET_READ_LENGTH           } from '../../modules/local/get_read_length'
include { GET_UNMAPPED              } from '../../modules/local/get_unmapped'
include { SPLIT_FASTQS              } from '../../modules/local/split_fastqs'
include { GET_ANCHORS               } from '../../modules/local/get_anchors'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { ANCHOR_SAMPLE_SCORES      } from '../../modules/local/anchor_sample_scores'
include { NORM_SCORES               } from '../../modules/local/norm_scores'

include { TRIMGALORE                } from '../../modules/nf-core/modules/trimgalore/main'

workflow ANALYZE_FASTQS {

    main:

    // Step 0: parse samplesheet
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

    // Step 0: definitions
    num_lines = params.chunk_size * params.n_iterations * 4
    num_chunk_lines = params.chunk_size * 4

    // Step 0: define lookahead parameter
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

    // Step 1: get significant anchors
    if (params.anchors_file) {
        // use input anchors file
        ch_anchors = params.anchors_file

    } else {

        /*
        // Trim fastqs
        */
        TRIMGALORE(
            ch_fastqs
        )

        /*
        // only use unmapped reads
        */

        if (params.unmapped) {

            GET_UNMAPPED(
                TRIMGALORE.out.fastq,
                params.index_bowtie
            )

            /*
            // Process to split each fastq into size read_chunk and gzip compress
            */

            SPLIT_FASTQS(
                GET_UNMAPPED.out.fastq,
                num_lines,
                num_chunk_lines
            )

            ch_split_fastqs = SPLIT_FASTQS.out.fastq.collect()
            }

        else {

            /*
            // Process to split each fastq into size read_chunk and gzip compress
            */

            SPLIT_FASTQS(
                TRIMGALORE.out.fastq,
                num_lines,
                num_chunk_lines,
                params.n_iterations
            )

            ch_split_fastqs = SPLIT_FASTQS.out.fastq.collect()

        }


        /*
        // Process to get list of candidate anchors via onthefly
        */

        GET_ANCHORS(
            ch_split_fastqs,
            params.n_iterations,
            params.chunk_size,
            params.kmer_size,
            file(params.input),
            params.target_counts_threshold,
            params.anchor_counts_threshold,
            params.anchor_freeze_threshold,
            params.read_freeze_threshold,
            params.anchor_score_threshold,
            params.anchor_mode,
            params.c_type,
            params.window_slide,
            lookahead,
            params.num_keep_anchors,
            params.use_std,
            params.compute_target_distance,
            params.max_ignorelist,
            params.distance_type,
            params.score_type
        )

        ch_anchors = GET_ANCHORS.out.anchors
    }

    // Step 2
    /*
    // Process to get consensus sequences and target counts for annchors
    */
    PARSE_ANCHORS(
        ch_fastqs,
        ch_anchors,
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

    // Step 3
    /*
    // Process to get anchor scores and anchor-target counts
    */
    ANCHOR_SAMPLE_SCORES(
        targets_samplesheet,
        params.bound_distance,
        params.max_distance,
        params.kmer_size,
        params.distance_type
    )

    anchor_sample_scores = ANCHOR_SAMPLE_SCORES.out.anchor_sample_scores.collect()
    anchor_target_counts = ANCHOR_SAMPLE_SCORES.out.anchor_target_counts.collect()

    // Step 4
    /*
    // Process to compute norm scores
    */
    NORM_SCORES(
        anchor_sample_scores,
        anchor_target_counts,
        file(params.input),
        params.use_std,
        params.kmer_size
    )

    emit:
    anchor_target_counts = anchor_target_counts

}
