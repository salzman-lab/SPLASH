

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

    // Process samplesheet
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

    // define number of lines to grab from fastq
    num_lines = params.chunk_size * params.n_iterations * 4
    num_chunk_lines = params.chunk_size * 4

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
    println(params.target_threshold)
    GET_ANCHORS(
        ch_split_fastqs,
        params.chunk_size,
        file(samplesheet),
        params.target_threshold,
        params.anchor_counts_threshold,
        params.anchor_freeze_threshold,
        params.anchor_mode,
        params.window_slide,
        params.read_len
    )

    ch_anchors = GET_ANCHORS.out.anchors

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
        params.adj_distance
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


