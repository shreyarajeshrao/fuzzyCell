d <- miniPipelineData()

test_that("confidenceWeights requires a SummarizedExperiment", {
    for (bad in list(getCellData(d$spe), NULL, "spe", 1)) {
        expect_error(confidenceWeights(bad, d$tree, d$conf),
            "`spe` must be a SummarizedExperiment, for example a SpatialExperiment")
    }
})

test_that("confidenceWeights reports cell types missing from the tree", {
    spe <- d$spe
    spe$cellType[1] <- "Ghost"
    expect_error(confidenceWeights(spe, d$tree, d$conf),
        "1 cell type \\(1 cell\\) in the data is not in cellTypeTree: 'Ghost'")
})

test_that("weighted rows hold the confidences, unweighted rows 0/1 labels", {
    w <- confidenceWeights(d$spe, d$tree, d$conf)
    u <- confidenceWeights(d$spe, d$tree, d$conf, weighted = FALSE)
    expect_equal(dim(w), dim(u))
    expect_true(all(u %in% c(0, 1)))
    inG <- d$spe$cellType %in% c("A", "B")
    conf <- d$conf$G[match(d$spe$cellID[inG], d$conf$G$cellID), ]
    expect_equal(unname(w[inG, "A"]), conf$A)
    expect_true(all(w[inG, "G"] == 1))
})
