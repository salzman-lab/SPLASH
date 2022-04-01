include { GET_FASTA             } from '../../modules/local/get_fasta'
include { BOWTIE2_ANNOTATION    } from '../../modules/local/bowtie2_annotation'
include { ADD_ANNOTATIONS       } from '../../modules/local/add_annotations'


workflow REANNOTATE {
    take:
    anchor_target_counts
    anchor_scores

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

    /*
    // Process to align anchors to each bowtie2 index
    */
    BOWTIE2_ANNOTATION(
        GET_FASTA.out.anchor_fasta,
        GET_FASTA.out.target_fasta,
        ch_indices
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
    // Process to merge anchor scores with anchor hits
    */
    ADD_ANNOTATIONS(
        anchor_hits_samplesheet,
        target_hits_samplesheet,
        anchor_scores
    )

}
