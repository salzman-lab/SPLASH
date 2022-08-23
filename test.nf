Channel
    .from(
        [1, ["A", "a"]],
        [1, ["B", "b"]],
        [-1, ["C", "c"]],
        [-1, ["D", "d"]]
    )

cj = Channel.from(1, -1)
samples = Channel
    .from(
        ["A", "a"],
        ["B", "b"],
        ["C", "c"],
        ["D", "d"]
    )

cj
    .combine(samples)
    .map{ it ->
        tuple(
            [it[1], it[2]],
            it[0]
        )
    }
    .groupTuple()
    .map{ it ->
        tuple(
            it[0],
            it[1].flatten()
        )
    }
    .view()
    .transpose()
    .view()
