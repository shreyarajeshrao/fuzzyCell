test_that("getCellData adds coordinates from a SpatialExperiment", {
    cells <- getCellData(miniSpe())
    expect_equal(nrow(cells), 4)
    expect_equal(cells$x, c(1, 2, 3, 4))
    expect_equal(cells$y, c(4, 3, 2, 1))
    expect_true(all(c("cellID", "imageID", "cellType") %in% names(cells)))
})

test_that("generateComparisons returns all pair types", {
    confList <- list("T cells" = NULL, "B lineage" = NULL)
    comp <- fuzzyCell:::generateComparisons(heartCellTypeTree, confList)
    expect_named(comp, c("groupGroup", "groupLeaf", "leafLeaf"))
    nLeaves <- nrow(buildCellTypeMap(heartCellTypeTree)$cellMap)
    expect_equal(nrow(comp$leafLeaf), nLeaves^2)
    expect_equal(nrow(comp$groupGroup), 4)
})

test_that("parallel random streams follow the session's seed", {
    draw <- function() BiocParallel::bplapply(1:3, function(i) stats::runif(1),
        BPPARAM = fuzzyCell:::parallelParam(1))
    set.seed(1)
    a <- draw()
    set.seed(1)
    expect_identical(draw(), a)
    expect_false(identical(draw(), a))
})

test_that("mixed models record singular fits without printing lmerTest's note",
    {
    skip_if_not_installed("lmerTest")
    # six subjects with two images each; the subjects explain none of the
    # variation, so the between-subject variance is estimated as zero
    dat <- data.frame(y = c(-0.63, 0.18, -0.84, 1.6, 0.33, -0.82, 2.49, 2.74,
        2.58, 1.69, 3.51, 2.39), outcome = rep(c("a", "b"), each = 6),
        .subject = rep(paste0("s", 1:6), each = 2))
    expect_no_message(fit <- fuzzyCell:::mixedModelTest(dat, "y", "outcome"))
    expect_equal(fit[["singular"]], 1)

    tab <- data.frame(cellType = c("A", "B", "C"), singular = c(1, 0, NA))
    expect_message(out <- fuzzyCell:::reportSingular(tab, "test proportions"),
        "for the test proportions: for 1 of the 2 models")
    expect_false("singular" %in% names(out))
    expect_no_message(fuzzyCell:::reportSingular(
        data.frame(singular = c(0, NA)), "x"))
})
