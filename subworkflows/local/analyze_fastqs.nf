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

    // use input anchors file
    ch_anchors = file(params.anchors_file)

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

}
