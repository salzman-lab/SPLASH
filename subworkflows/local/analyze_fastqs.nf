include { GET_READ_LENGTH           } from '../../modules/local/get_read_length'
include { GET_UNMAPPED              } from '../../modules/local/get_unmapped'
include { FETCH_ANCHORS             } from '../../modules/local/fetch_anchors'
include { COUNT_ANCHORS             } from '../../modules/local/count_anchors'
include { STRATIFY_ANCHORS          } from '../../modules/local/stratify_anchors'
include { GET_ABUNDANT_ANCHORS      } from '../../modules/local/get_abundant_anchors'
include { MERGE_ABUNDANT_ANCHORS    } from '../../modules/local/merge_abundant_anchors'
include { GET_ANCHORS_AND_SCORES    } from '../../modules/local/get_anchors_and_scores'
include { MERGE_ANCHOR_SCORES       } from '../../modules/local/merge_anchor_scores'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { MERGE_TARGET_COUNTS       } from '../../modules/local/merge_target_counts'

include { TRIMGALORE                } from '../../modules/nf-core/modules/trimgalore/main'

workflow ANALYZE_FASTQS {

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

    if (! params.skip_trimming) {
        /*
        // Trim fastqs
        */
        TRIMGALORE(
            ch_fastqs
        )
        ch_fastqs = TRIMGALORE.out.fastq
    }

    /*
    // Check if we are only using unmapped reads
    */
    if (params.unmapped) {
        /*
        // Get unmapped reads
        */
        GET_UNMAPPED(
            ch_fastqs,
            params.index_bowtie
        )

        ch_fastqs = GET_UNMAPPED.out.fastq

    }

    /*
    // Process to get all candidate anchors and targets
    */
    FETCH_ANCHORS(
        ch_fastqs,
        params.run_type,
        params.num_lines,
        params.kmer_size,
        lookahead,
        params.anchor_mode,
        params.window_slide

    )

    /*
    // Process to count anchors and targets
    */
    COUNT_ANCHORS(
        FETCH_ANCHORS.out.seqs
    )

    /*
    // Process to stratify counts files by 3mers
    */
    STRATIFY_ANCHORS(
        COUNT_ANCHORS.out.seqs.collect(),
        params.stratify_level
    )

    /*
    // Process to filter kmer counts for abundant anchors
    */
    GET_ABUNDANT_ANCHORS(
        STRATIFY_ANCHORS.out.seqs.flatten(),
        params.anchor_count_threshold
    )

    /*
    // Process to get significant anchors and their scores
    */
    GET_ANCHORS_AND_SCORES(
        GET_ABUNDANT_ANCHORS.out.seqs,
        file(params.input),
        params.distance_type,
        params.max_targets,
        params.max_dist,
        params.bonfer,
        params.pval_threshold,
        params.run_type
    )

    // Merge all scores from all slices and output
    GET_ANCHORS_AND_SCORES.out.scores
        .collectFile(
            name: 'scores.tsv',
            storeDir: "${params.outdir}",
            keepHeader: true,
            skip: 1
        )
        .set{anchor_scores}

    // Merge all anchors from all slices and output
    GET_ANCHORS_AND_SCORES.out.anchors
        .collectFile(
            name: 'all_anchors.tsv',
            keepHeader: true,
            skip: 1
        )
        .set{all_anchors}

    /*
    // Process to output top 5000 anchors as sorted by pvalue
    */
    MERGE_ANCHOR_SCORES(
        all_anchors
    )

    /*
    // Process to get consensus sequences and target counts for annchors
    */
    PARSE_ANCHORS(
        ch_fastqs,
        MERGE_ANCHOR_SCORES.out.anchors.first(),
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

    emit:
    anchor_target_counts    = MERGE_TARGET_COUNTS.out.anchor_target_counts.first()
    anchor_scores           = anchor_scores
    ch_consensus_fasta      = PARSE_ANCHORS.out.consensus_fasta.collect()

}
