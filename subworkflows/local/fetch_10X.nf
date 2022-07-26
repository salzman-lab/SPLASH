include { UMITOOLS              } from '../../modules/local/umitools'
include { EXTRACT_CBC_UMI       } from '../../modules/local/extract_cbc_umi'
include { DEDUP_CBC_UMI         } from '../../modules/local/dedup_cbc_umi'

workflow FETCH_10X {

    take:
    ch_paired_fastqs

    main:

    UMITOOLS(
        ch_paired_fastqs
    )

    EXTRACT_CBC_UMI(
        UMITOOLS.out.fastqs,
        params.num_reads_first_pass
    )

    DEDUP_CBC_UMI(
        EXTRACT_CBC_UMI.out.seqs
    )

    emit:
    fastqs = DEDUP_CBC_UMI.out.seqs

}
