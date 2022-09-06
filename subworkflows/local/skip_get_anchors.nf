include { GET_READ_LENGTH           } from '../../modules/local/get_read_length'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { MERGE_TARGET_COUNTS       } from '../../modules/local/merge_target_counts'
include { GET_FASTA                 } from '../../modules/local/get_fasta'
include { BOWTIE2_ANNOTATION        } from '../../modules/local/bowtie2_annotation'
include { MERGE_ANNOTATIONS         } from '../../modules/local/merge_annotations'
include { GENOME_ALIGNMENT      } from '../../modules/local/genome_alignment'
include { GENOME_ANNOTATIONS    } from '../../modules/local/genome_annotations'
include { STAR_ALIGN            } from '../../modules/local/star_align'

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

    /*
    // Process to get consensus sequences and target counts for annchors
    */
    PARSE_ANCHORS(
        ch_fastqs,
        file(params.anchors_file),
        params.num_parse_anchors_reads,
        params.consensus_length,
        params.kmer_size,
        params.target_size,
        params.direction,
        lookahead
    )

    // Create samplesheet of target counts files
    PARSE_ANCHORS.out.targets
        .collectFile() { file ->
            def X=file; X.toString() + '\n'
        }
        .set{ targets_samplesheet }

    ch_consensus_fastas = PARSE_ANCHORS.out.consensus_fasta

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

    ch_anchor_target_fastas = GET_FASTA.out.fasta

    /*
    // Cartesian join of anchor+target fastas and all bowtie2 indices
    */
    ch_anchor_target_indices = ch_anchor_target_fastas
        .flatten()
        .combine(ch_indices)

    /*
    // Process to align anchors to each bowtie2 index
    */
    BOWTIE2_ANNOTATION(
        ch_anchor_target_indices
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

    /*
    // Make one channel containing anchor, target, and consensus fastas
    */
    ch_anchor_target_fastas
        .flatten()
        .mix(ch_consensus_fastas)
        .flatten()
        .set{ch_fastas}

    /*
    // Process to align targets and anchors to genome
    */
    GENOME_ALIGNMENT(
        ch_fastas,
        params.genome_index,
        params.transcriptome_index
    )

    /*
    // Process to run gene and exon annotations
    */
    GENOME_ANNOTATIONS(
        GENOME_ALIGNMENT.out.bam_tuple,
        params.gene_bed,
        params.exon_starts_bed,
        params.exon_ends_bed
    )

}
