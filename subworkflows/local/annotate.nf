include { ELEMENT_ALIGNMENT         } from '../../modules/local/element_alignment'
include { ELEMENT_ANNOTATIONS       } from '../../modules/local/element_annotations'
include { SUMMARIZE                 } from '../../modules/local/summarize'
include { GENOME_ALIGNMENT          } from '../../modules/local/genome_alignment'
include { GENOME_ANNOTATIONS        } from '../../modules/local/genome_annotations'
include { PREPARE_CONSENSUS         } from '../../modules/local/prepare_consensus'
include { MERGE_CONSENSUS           } from '../../modules/local/merge_consensus'
include { CONSENSUS_ALIGNMENT       } from '../../modules/local/consensus_alignment'
include { SPLICING_ANNOTATIONS      } from '../../modules/local/splicing_annotations'
include { ADDITIONAL_SUMMARY        } from '../../modules/local/additional_summary'


workflow ANNOTATE {

    take:
    anchor_scores
    anchor_target_counts
    ch_consensus_fastas
    ch_anchor_target_fastas

    main:

    // create samplesheet of bowtie2 indices
    ch_indices = Channel.fromPath(params.element_annotations_samplesheet)
        .splitCsv(
            header: false
        )
        .map{ row ->
            row[0]
        }


    /*
    // Cartesian join of anchor+target fastas and all bowtie2 indices
    */
    ch_anchor_target_indices = ch_anchor_target_fastas
        .flatten()
        .combine(ch_indices)

    /*
    // Process to align anchors to each bowtie2 index
    */
    ELEMENT_ALIGNMENT(
        ch_anchor_target_indices
    )

    /*
    // Process to merge scores with hits
    */
    ELEMENT_ANNOTATIONS(
        ELEMENT_ALIGNMENT.out.hits.collect()
    )

    /*
    // Process to run postprocessing annotations
    */
    SUMMARIZE(
        anchor_scores,
        anchor_target_counts,
        ELEMENT_ANNOTATIONS.out.annotated_anchors,
        ELEMENT_ANNOTATIONS.out.annotated_targets,
        params.run_blast
    )

    /*
    // Make one channel containing anchor, target, and consensus fastas
    */
    ch_fastas = ch_anchor_target_fastas
        .flatten()
        .mix(ch_consensus_fastas)
        .flatten()

    if (params.run_annotations) {
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
        // Process to prepare consensus fastas for one STAR alignment
        */
        PREPARE_CONSENSUS(
            ch_consensus_fastas.flatten()
        )

        /*
        // Process to concatenate consensus fastas for one STAR alignment
        */
        MERGE_CONSENSUS(
            PREPARE_CONSENSUS.out.fasta.collect()
        )

        fasta = MERGE_CONSENSUS.out.fasta
        /*
        // Process to get splice junctions with STAR
        */
        CONSENSUS_ALIGNMENT(
            fasta,
            params.star_index,
            params.gtf
        )

        GENOME_ANNOTATIONS.out.annotations
            .filter{
                file ->
                file.name.contains('genome_annotations_anchor.tsv')
            }
            .set{genome_annotations_anchors}

        /*
        // Process to get called exons from bam file
        */
        SPLICING_ANNOTATIONS(
            CONSENSUS_ALIGNMENT.out.bam,
            CONSENSUS_ALIGNMENT.out.unmapped_fasta,
            params.gene_bed,
            params.ann_AS_gtf,
            fasta,
            genome_annotations_anchors
        )

        /*
        // Process to make additional summary file
        */
        ADDITIONAL_SUMMARY(
            SPLICING_ANNOTATIONS.out.consensus_called_exons,
            SUMMARIZE.out.tsv
        )
    }

}
