include { GET_READ_LENGTH           } from '../../modules/local/get_read_length'
include { SPLIT_FASTQS              } from '../../modules/local/split_fastqs'
include { GET_ANCHORS               } from '../../modules/local/get_anchors'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { COMPUTE_ANCHOR_SCORES     } from '../../modules/local/compute_anchor_scores'
include { NORM_SCORES               } from '../../modules/local/norm_scores'

include { TRIMGALORE                } from '../../modules/nf-core/modules/trimgalore/main'

workflow ANALYZE_FASTQS {

    main:

    // parse samplesheet
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

    // definitions
    num_lines = params.chunk_size * params.n_iterations * 4
    num_chunk_lines = params.chunk_size * 4

    if (params.use_read_length) {
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
        // Trim fastqs
        */
        TRIMGALORE(
            ch_fastqs
        )

        /*
        // Process to split each fastq into size read_chunk and gzip compress
        */
        SPLIT_FASTQS(
            TRIMGALORE.out.fastq,
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
            file(params.input),
            params.target_counts_threshold,
            params.anchor_counts_threshold,
            params.anchor_freeze_threshold,
            params.read_freeze_threshold,
            params.anchor_score_threshold,
            params.anchor_mode,
            params.c_type,
            params.window_slide,
            looklength,
            params.num_keep_anchors,
            params.use_std,
            params.compute_target_distance,
            params.max_ignorelist
        )

        ch_anchors = GET_ANCHORS.out.anchors
    }

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
        targets_samplesheet,
        params.bound_distance,
        params.max_distance,
        params.kmer_size
    )

    anchor_scores           = COMPUTE_ANCHOR_SCORES.out.anchor_scores.first()
    anchor_target_counts    = COMPUTE_ANCHOR_SCORES.out.anchor_target_counts.first()

    /*
    // Process to compute norm scores
    */
    NORM_SCORES(
        anchor_scores,
        anchor_target_counts,
        file(params.input),
        params.kmer_size
    )

    emit:
    anchor_target_counts = anchor_target_counts

}
