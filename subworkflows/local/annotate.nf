include { GET_FASTA             } from '../../modules/local/get_fasta'
include { BOWTIE2_ANNOTATION    } from '../../modules/local/bowtie2_annotation'
include { MERGE_ANNOTATIONS     } from '../../modules/local/merge_annotations'
include { POSTPROCESSING        } from '../../modules/local/postprocessing'
include { GENOME_ALIGNMENT      } from '../../modules/local/genome_alignment'
include { GENOME_ANNOTATIONS    } from '../../modules/local/genome_annotations'
include { STAR_ALIGN            } from '../../modules/local/star_align'
include { ANNOTATE_SPLICES      } from '../../modules/local/annotate_splices'


workflow ANNOTATE {
    take:
    anchor_target_counts
    anchor_scores
    ch_consensus_fastas

    main:

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
        anchor_target_counts
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
    // Process to run postprocessing annotations
    */
    POSTPROCESSING(
        anchor_scores,
        anchor_target_counts,
        MERGE_ANNOTATIONS.out.annotated_anchors,
        MERGE_ANNOTATIONS.out.annotated_targets,
        params.run_blast
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

    /*
    // Process to get splice junctions with STAR
    */
    STAR_ALIGN(
        ch_consensus_fastas.flatten(),
        params.star_index,
        params.gtf
    )

    // /*
    // // Process to annotate splice junctions
    // */
    // ANNOTATE_SPLICES(
    //     STAR_ALIGN.out.junctions,
    //     params.exon_pickle,
    //     params.splice_pickle
    // )
}
