include { UMITOOLS              } from '../../modules/local/umitools'
include { EXTRACT_CBC_UMI       } from '../../modules/local/extract_cbc_umi'
include { DEDUP_CBC_UMI         } from '../../modules/local/dedup_cbc_umi'
include { CONCAT_FASTQS         } from '../../modules/local/concat_fastqs'

workflow PREPROCESS_10X {

    take:
    ch_paired_fastqs

    main:

    // Run umitools first if that has not been done already
    if (params.run_umitools) {

        /*
        // Process to run umitools whitelist and extract
        */
        UMITOOLS(
            ch_paired_fastqs
        )
        ch_umitools_fastqs = UMITOOLS.out.fastqs

    } else {
        ch_umitools_fastqs = ch_paired_fastqs
            .map{ it ->
                tuple(
                    it[0],
                    it[2]
                )
            }

    }

    if (params.cells_split_across_lanes) {
        // We use this case if the same CBC are split across multiple lanes

        // Group files by channel, for deduplication by channel
        ch_umitools_fastqs
            .map { it ->
                tuple(
                    it[0],
                    [it[1]]
                )
            }
            .groupTuple()
            .map { it ->
                tuple(
                    it[0],
                    it[1].flatten()
                )
            }
            .view()
            .set { ch_merged_fastqs }

        /*
        // Process to concat fastqs by sample if cells are split across lanes
        */
        CONCAT_FASTQS(
            ch_merged_fastqs
        )

        ch_channel_fastqs = CONCAT_FASTQS.out.fastqs

    } else {
        ch_channel_fastqs = ch_umitools_fastqs
    }

    /*
    // Process to extract CBC and UMI from fastq header, and write out CBC_UMI with read
    */
    EXTRACT_CBC_UMI(
        ch_channel_fastqs,
        params.num_fastq_reads
    )

    /*
    // Process to deduplicate on sample_CBC_UMI
    */
    DEDUP_CBC_UMI(
        EXTRACT_CBC_UMI.out.seqs
    )

    emit:
    fastqs = DEDUP_CBC_UMI.out.seqs

}
