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
    anchors_pvals
    anchor_target_counts
    ch_consensus_fastas
    ch_anchor_target_fastas

    main:

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

    /*
    // Cartesian join of anchor+target fastas and all bowtie2 indices
    */
    ch_anchor_target_indices = ch_anchor_target_fastas
        .flatten()
        .combine(ch_indices)

    /*
    // Process: Align anchors to each bowtie2 index
    */
    ELEMENT_ALIGNMENT(
        ch_anchor_target_indices
    )

    /*
    // Process: Merge scores with hits
    */
    ELEMENT_ANNOTATIONS(
        ELEMENT_ALIGNMENT.out.hits.collect()
    )

    /*
    // Process: Run postprocessing annotations
    */
    SUMMARIZE(
        anchors_pvals,
        anchor_target_counts,
        ELEMENT_ANNOTATIONS.out.annotated_anchors,
        ELEMENT_ANNOTATIONS.out.annotated_targets
    )

    if (params.run_annotations) {
        /*
        // Process: Align targets and anchors to genome
        */
        GENOME_ALIGNMENT(
            ch_anchor_target_fastas.flatten(),
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

        /*
        // Process: Prepare consensus fastas for one STAR alignment
        */
        PREPARE_CONSENSUS(
            ch_consensus_fastas.flatten()
        )

        /*
        // Process: Concatenate consensus fastas for one STAR alignment
        */
        MERGE_CONSENSUS(
            PREPARE_CONSENSUS.out.fasta.collect()
        )

        consensus_fasta = MERGE_CONSENSUS.out.fasta
        /*
        // Process: Get splice junctions with STAR
        */
        CONSENSUS_ALIGNMENT(
            consensus_fasta,
            params.star_index,
            params.gtf
        )

        genome_annotations_anchors = GENOME_ANNOTATIONS.out.annotated_anchors

        /*
        // Process: Get called exons from bam file
        */
        SPLICING_ANNOTATIONS(
            CONSENSUS_ALIGNMENT.out.bam,
            CONSENSUS_ALIGNMENT.out.unmapped_fasta,
            params.gene_bed,
            consensus_fasta,
            genome_annotations_anchors
        )

        /*
        // Process: Make additional summary file
        */
        ADDITIONAL_SUMMARY(
            SPLICING_ANNOTATIONS.out.consensus_called_exons,
            SUMMARIZE.out.tsv
        )

        additional_summary = ADDITIONAL_SUMMARY.out.tsv


    } else {
        additional_summary          = null
        genome_annotations_anchors  = null

    }

    emit:
    additional_summary          = additional_summary
    genome_annotations_anchors  = genome_annotations_anchors
}
