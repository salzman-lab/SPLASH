include { MAKE_FASTA                } from '../../modules/local/make_fasta'
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

    main:

    // Make fasta from anchors file
    MAKE_FASTA(
        file(params.anchors)
    )

    // Set anchors channel
    ch_anchors = MAKE_FASTA.out.fasta

    // Parse samplesheet of bowtie2 indices
    element_annotations_samplesheet = file(params.element_annotations_samplesheet, checkIfExists: true)
    ch_indices = Channel.fromPath(element_annotations_samplesheet)
        .splitCsv(
            header: false
        )
        .map{ row ->
            row[0]
        }

    /*
    // Cartesian join of anchors and all bowtie2 indices
    */
    ch_anchor_indices = ch_anchors
        .flatten()
        .combine(ch_indices)

    /*
    // Process to align anchors to each bowtie2 index
    */
    ELEMENT_ALIGNMENT(
        ch_anchor_indices
    )

    /*
    // Process to merge scores with hits
    */
    ELEMENT_ANNOTATIONS(
        ELEMENT_ALIGNMENT.out.hits.collect()
    )

    /*
    // Process to align targets and anchors to genome
    */
    GENOME_ALIGNMENT(
        ch_anchors,
        params.genome_index,
        params.transcriptome_index
    )

    /*
    // Process to run gene and exon annotations
    */
    GENOME_ANNOTATIONS(
        GENOME_ALIGNMENT.out.bam_tuple,
        params.gene_bed
    )

    /*
    // Process to run postprocessing annotations
    */
    SUMMARIZE(
        file(params.anchors),
        ELEMENT_ANNOTATIONS.out.anchors,
        GENOME_ANNOTATIONS.out.anchors
    )

}
