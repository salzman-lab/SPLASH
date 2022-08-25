ch_fastqs = Channel
    .fromPath(file("samplesheet_test.csv"))
    .splitCsv(
        header: false
    )
    .map { row ->
        tuple(
            row[1],
            file(row[2])
        )
    }

l = ch_fastqs
    .map{ it[1]}
    .toSortedList()
    .view()

Channel.fromList(l).view()
// Channel
//     .fromList(
//         ch_fastqs
//             .map{ it[1]}
//             .toSortedList()
//     )
//     .view()
