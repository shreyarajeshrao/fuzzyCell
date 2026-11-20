# k-means clusters on LISA curves; features that are NA (a `to` type absent
# from the image) are set to 0, i.e. no attraction or repulsion
clusterNiches <- function(lisaMat, k, ...) {
    lisaMat[is.na(lisaMat)] <- 0
    km <- kmeans(lisaMat, centers = k, ...)
    setNames(paste0("niche", km$cluster), rownames(lisaMat))
}

#' Niche identification
#'
#' Defines niches by k-means clustering of per-cell feature vectors, once with
#' hard cell type membership and once with bootstrap-confidence weights. The
#' features of each cell are its local indicators of spatial association from
#' lisaClust (`lisaClust::lisa()`) at several radii.
#'
#' @details
#' For every cell, every label in `toLabels` and every radius in `r`, the
#' feature is the local co-localisation of the cell's neighbours with that label
#' (lisaClust's local K or L function, centred so that 0 means no attraction or
#' repulsion). The weighted features count each neighbour with its bootstrap
#' confidence of having the label, instead of 0 or 1 (see
#' [confidenceWeights()]), and a label's density is
#' its summed weight per unit area.
#'
#' Niche names (`niche1`, `niche2`, ...) are arbitrary. K-means starts
#' from random cluster centres, so the niches depend on the random seed; call
#' [set.seed()] first for reproducible niches. The
#' unweighted and weighted clusterings both start from the same random state, so
#' that if the weights made no difference the two sets of niches would be
#' identical.
#'
#' @inheritParams runColocalisationModels
#' @param toLabels Labels (leaves and/or groups) to build features from;
#'   defaults to every leaf and group in the tree. An error is raised for a
#'   label that is not in the tree.
#' @param r Radii of the local curves: one or more positive numbers.
#' @param k Number of niches, a whole number of at least 2 and below the number
#'   of cells.
#' @param window The area of each image that densities are computed over:
#'   `"convex"` (the default), the convex hull of the image's cells, or
#'   `"square"`, the rectangle around them, as in
#'   [runColocalisationModels()].
#' @param lisaFunc The local function, `"K"` or `"L"`, as in
#'   `lisaClust::lisa()`.
#'
#' @return A list with:
#' * `lisaMatUnweighted`, `lisaMatWeighted`: the features, matrices
#'   with one row per cell (named by cell ID) and one column per label in
#'   `toLabels` and radius in `r`;
#' * `nicheUnweighted`, `nicheWeighted`: each cell's niche
#'   (`"niche1"`, `"niche2"`, ...), as character vectors named by cell ID.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' set.seed(1)
#' niches <- runNicheIdentification(lungSubset, lungCellTypeTree, lungConfList,
#'     r = 50, k = 6)
#' table(weighted = niches$nicheWeighted, unweighted = niches$nicheUnweighted)
#' @name nicheIdentification
#' @export
runNicheIdentification <- function(spe, cellTypeTree, confList,
        toLabels = NULL, r = c(20, 50, 100), k, window = c("convex", "square"),
        lisaFunc = "K", cores = 1, cellType = "cellType", cellID = "cellID",
        imageID = "imageID", x = "x", y = "y") {
    # check the arguments before doing any work
    window <- match.arg(window)
    checkSpe(spe)
    if (missing(k)) {
        stop("`k`, the number of niches, must be given", call. = FALSE)
    }
    checkNumberOfClusters(k)
    checkRadii(r)
    toLabels <- nicheLabels(toLabels, cellTypeTree)

    cells <- standardiseCells(getCellData(spe), list(cellID = cellID,
        cellType = cellType, imageID = imageID, x = x, y = y))
    if (k >= nrow(cells)) {
        stop("`k` = ", k, " must be smaller than the number of cells (",
            nrow(cells), ")", call. = FALSE)
    }

    lisaMat <- lapply(c(unweighted = FALSE, weighted = TRUE), function(wtd) {
        message("Computing ", if (wtd) "weighted" else "unweighted",
            " LISA curves...")
        w <- weightMatrix(cells, cellTypeTree, confList, wtd)
        lisaClust::lisa(cells[, c("cellID", "imageID", "cellType", "x", "y")],
            r = r, window = window, lisaFunc = lisaFunc, cores = cores,
            cellWeights = w[, intersect(colnames(w), toLabels), drop = FALSE])
    })

    # both clusterings start from the same random state, so that they differ
    # only through the features and not through different random starts. The
    # draw makes sure the session has a random state to keep.
    message("Clustering niches...")
    stats::runif(1)
    nicheUnweighted <- withr::with_preserve_seed(
        clusterNiches(lisaMat$unweighted, k = k))
    nicheWeighted <- clusterNiches(lisaMat$weighted, k = k)

    list(
        lisaMatUnweighted = lisaMat$unweighted,
        lisaMatWeighted = lisaMat$weighted,
        nicheUnweighted = nicheUnweighted, nicheWeighted = nicheWeighted
    )
}

# The labels to build niche features from: `toLabels`, checked against the
# tree, or by default every cell type and group of the tree.
nicheLabels <- function(toLabels, cellTypeTree) {
    if (!is.null(toLabels)) {
        checkTreeLabels(toLabels, cellTypeTree)
        return(toLabels)
    }
    treeInfo <- buildCellTypeMap(cellTypeTree)
    union(treeInfo$cellMap$cellType, treeInfo$groupsWithChildren)
}

#' Test niche proportions against an outcome
#'
#' Computes the proportion of each image's cells in each niche and tests every
#' niche against a two-level outcome with a Wilcoxon or t-test.
#'
#' @details
#' The `estimate` is the difference in the niche proportion for the second
#' level of the outcome compared with the first (alphabetical order for text,
#' factor order for factors), as in
#' [runColocalisationModels()]. For
#' `test = "wilcox"` this is the Hodges-Lehmann estimate of the shift in
#' location.
#'
#' Each image is one observation. When `subject` is given and some subject has
#' several images, those images are not independent:
#' * with `test = "ttest"`, each niche is tested with a linear mixed model
#'   with a random intercept per subject (`lmerTest::lmer()`, as spicyR fits for
#'   co-localisation), which needs the lmerTest package;
#' * with `test = "wilcox"`, each subject's cells are pooled into one set
#'   of niche proportions, and the subjects are compared (as `spicyR::colTest()`
#'   does). This needs every image of a subject to have the same outcome.
#'
#' A message says which is used.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   a cell ID, an image and the outcome for every cell in `colData`.
#' @param cellID,imageID Names of the `colData(spe)` columns holding a
#'   unique cell ID (matching the names of the niche vectors) and the image
#'   each cell is in.
#' @param nicheUnweighted,nicheWeighted Named cell-to-niche vectors, as
#'   returned by [runNicheIdentification()].
#' @param outcomeCol Name of the per-image outcome column in `colData(spe)`.
#'   It must have exactly two levels (subset your data to two conditions first)
#'   and a single value per image. An error is raised otherwise, before any
#'   computation.
#' @param test `"wilcox"` or `"ttest"`.
#' @param subject Name of the `colData(spe)` column identifying the subject
#'   (for example the patient) each image comes from, or `NULL`. If given,
#'   every image must belong to one subject; subjects with several images are
#'   handled as described in Details. With `NULL` every image is treated as a
#'   separate subject, and a message says so.
#'
#' @return A list with `testUnweighted`, `testWeighted` (one row per niche) and
#'   `propUnweighted`, `propWeighted` (per image niche proportions, to pass
#'   to [nicheProportionBoxplot()]).
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' set.seed(1)
#' niches <- runNicheIdentification(lungSubset, lungCellTypeTree, lungConfList,
#'     r = 50, k = 6)
#' res <- testNicheOutcome(lungSubset, niches$nicheUnweighted,
#'     niches$nicheWeighted, outcomeCol = "disease_status", subject = "subject")
#' res$testWeighted
#' @export
testNicheOutcome <- function(spe, nicheUnweighted, nicheWeighted, outcomeCol,
        test = c("wilcox", "ttest"), subject = NULL, cellID = "cellID",
        imageID = "imageID") {
    test <- match.arg(test)

    # check the arguments before doing any work
    checkSpe(spe)
    cells <- standardiseCells(getCellData(spe), list(cellID = cellID,
        imageID = imageID))
    checkNicheVector(nicheUnweighted, "nicheUnweighted", cells)
    checkNicheVector(nicheWeighted, "nicheWeighted", cells)
    checkOutcome(cells, outcomeCol, subject)

    cells$nicheUnweighted <- nicheUnweighted[match(cells$cellID,
        names(nicheUnweighted))]
    cells$nicheWeighted <- nicheWeighted[match(cells$cellID,
        names(nicheWeighted))]
    propUnweighted <- nicheProportions(cells, "nicheUnweighted")
    propWeighted <- nicheProportions(cells, "nicheWeighted")
    # recorded so that nicheProportionBoxplot() can say which table it was given
    attr(propUnweighted, "nicheType") <- "unweighted"
    attr(propWeighted, "nicheType") <- "weighted"

    # the unit each test is run on: images, or subjects whose images are pooled
    units <- nicheTestUnits(cells, outcomeCol, subject, test)
    testProps <- function(niche, perImage) {
        if (units$pooled) nicheProportions(units$cells, niche) else perImage
    }
    list(
        testUnweighted = nicheTestTable(testProps("nicheUnweighted",
            propUnweighted), units$outcomeDat, outcomeCol, units$test,
            "unweighted niche proportions"),
        testWeighted = nicheTestTable(testProps("nicheWeighted",
            propWeighted), units$outcomeDat, outcomeCol, units$test,
            "weighted niche proportions"),
        propUnweighted = propUnweighted, propWeighted = propWeighted
    )
}

# How the niche proportions are tested. With one image per subject (or no
# `subject`), each image is one observation for `test`. When some subject has
# several images, "ttest" becomes a linear mixed model with a random intercept
# per subject ("mixed"), and "wilcox" pools each subject's cells into one
# observation. Returns the test to run, whether images are pooled, the cells to
# compute proportions from, and the outcome of each unit (`outcomeDat`, with
# the unit's ID in the imageID column and its subject in `.subject`).
nicheTestUnits <- function(cells, outcomeCol, subject, test) {
    keep <- c("imageID", outcomeCol)
    if (!is.null(subject)) {
        cells$.subject <- as.character(cells[[subject]])
        keep <- c(keep, ".subject")
    }
    perImage <- unique(cells[, keep])
    if (!hasRepeatedImages(cells, subject)) {
        return(list(test = test, pooled = FALSE, cells = cells,
            outcomeDat = perImage))
    }
    if (test == "ttest") {
        needLmerTest()
        message("Some subjects have several images: fitting a linear mixed ",
            "model with a random intercept per subject in place of the ",
            "t-test.")
        return(list(test = "mixed", pooled = FALSE, cells = cells,
            outcomeDat = perImage))
    }
    outcomes <- tapply(perImage[[outcomeCol]], perImage$.subject,
        function(v) length(unique(v)))
    if (any(outcomes > 1)) {
        stop("Some subjects have images with different outcomes (e.g. '",
            names(outcomes)[outcomes > 1][1], "'), so their images cannot be ",
            "pooled for the Wilcoxon test. Use test = \"ttest\", which fits a ",
            "mixed model.", call. = FALSE)
    }
    message("Some subjects have several images: pooling each subject's cells ",
        "into one set of niche proportions for the Wilcoxon test.")
    pooledCells <- cells
    pooledCells$imageID <- pooledCells$.subject
    perSubject <- unique(pooledCells[, c("imageID", outcomeCol)])
    list(test = test, pooled = TRUE, cells = pooledCells,
        outcomeDat = perSubject)
}

# Checks that `vec`, the argument called `arg`, is a vector of niches named by
# cell ID, with at least some of the cell IDs of `cells`.
checkNicheVector <- function(vec, arg, cells) {
    if (is.null(names(vec)) || !any(cells$cellID %in% names(vec))) {
        stop("`", arg, "` must be a vector of niches named by cell ID, with ",
            "at least some of the cell IDs in `cells`", call. = FALSE)
    }
    invisible(TRUE)
}

# The proportion of each image's cells in each niche, in long format (imageID,
# niche, proportion). `niche` names the column of `cells` holding the niches.
nicheProportions <- function(cells, niche) {
    props <- table(sample = cells$imageID, cluster = cells[[niche]])
    props <- props / rowSums(props)
    data.frame(
        imageID = rep(rownames(props), times = ncol(props)),
        niche = rep(colnames(props), each = nrow(props)),
        proportion = as.vector(props)
    )
}

# Tests every niche's proportion against the outcome with a Wilcoxon test, a
# t-test or a linear mixed model (`test` "wilcox", "ttest" or "mixed"), and
# returns one row per niche (estimate, statistic, p-value and BH-adjusted
# p-value), ordered by p-value. `propDat` is from nicheProportions(), and
# `outcomeDat` has one row per image (imageID, the outcome column and, for the
# mixed model, `.subject`).
nicheTestTable <- function(propDat, outcomeDat, outcomeCol, test,
        label = "niche proportions") {
    dat <- propDat
    rows <- match(dat$imageID, outcomeDat$imageID)
    dat[[outcomeCol]] <- outcomeDat[[outcomeCol]][rows]
    if (test == "mixed") dat$.subject <- outcomeDat$.subject[rows]
    formula <- as.formula(paste0("proportion ~ `", outcomeCol, "`"))
    noTest <- function(niche) {
        data.frame(niche = niche, estimate = NA, statistic = NA,
            p.value = NA, singular = NA)
    }

    rows <- lapply(unique(dat$niche), function(niche) {
        sub <- dat[dat$niche == niche & !is.na(dat$proportion) &
            !is.na(dat[[outcomeCol]]), , drop = FALSE]
        if (length(unique(sub[[outcomeCol]])) != 2 || nrow(sub) < 5) {
            return(noTest(niche))
        }
        if (test == "mixed") {
            fit <- mixedModelTest(sub, "proportion", outcomeCol)
            if (is.null(fit)) return(noTest(niche))
            return(data.frame(niche = niche, estimate = fit[["estimate"]],
                statistic = fit[["statistic"]], p.value = fit[["p.value"]],
                singular = fit[["singular"]]))
        }
        fit <- twoGroupTest(sub, formula, test)
        if (is.null(fit)) return(noTest(niche))
        data.frame(niche = niche, estimate = fit[["estimate"]],
            statistic = fit[["statistic"]], p.value = fit[["p.value"]],
            singular = NA)
    })
    out <- reportSingular(do.call(rbind, rows), label)
    out$p.adj <- p.adjust(out$p.value, method = "BH")
    out <- out[order(out$p.value), , drop = FALSE]
    rownames(out) <- NULL
    out
}

# A Wilcoxon test or a (Welch) t-test of the proportion on the outcome in
# `sub`: the estimate (the second level minus the first), the statistic and the
# p-value, or NULL if the test fails.
twoGroupTest <- function(sub, formula, test) {
    fit <- tryCatch({
        if (test == "wilcox") {
            stats::wilcox.test(formula, data = sub, conf.int = TRUE)
        } else {
            stats::t.test(formula, data = sub)
        }
    }, error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    # t.test() gives the two group means, so the difference is taken here;
    # wilcox.test() estimates first minus second, so its sign is flipped
    estimate <- if (test == "wilcox") {
        -unname(fit$estimate)
    } else {
        unname(diff(fit$estimate))
    }
    c(estimate = estimate, statistic = unname(fit$statistic),
        p.value = fit$p.value)
}

# Checks that `k`, the number of clusters, is a single whole number of at
# least 2.
checkNumberOfClusters <- function(k, arg = "k") {
    if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != round(k) ||
        k < 2) {
        stop("`", arg, "` (the number of niches) must be a single whole ",
            "number of at least 2", call. = FALSE)
    }
    invisible(k)
}

# Checks that every label is a cell type or a group of the tree.
checkTreeLabels <- function(labels, cellTypeTree, arg = "toLabels") {
    treeInfo <- buildCellTypeMap(cellTypeTree)
    valid <- union(treeInfo$cellMap$cellType, treeInfo$groupsWithChildren)
    if (!is.character(labels) || length(labels) == 0 || anyNA(labels)) {
        stop("`", arg, "` must be NULL or a character vector of cell types ",
            "and groups in cellTypeTree", call. = FALSE)
    }
    unknown <- setdiff(labels, valid)
    if (length(unknown) > 0) {
        shownUnknown <- paste0("'", utils::head(unknown, 5), "'",
            collapse = ", ")
        shownValid <- paste(utils::head(valid, 10), collapse = ", ")
        stop("`", arg, "` has ", length(unknown), " label",
            if (length(unknown) > 1) "s", " that ",
            if (length(unknown) > 1) "are" else "is",
            " not a cell type or group in cellTypeTree: ", shownUnknown,
            if (length(unknown) > 5) ", ...", ". Available labels: ",
            shownValid, if (length(valid) > 10) ", ...", call. = FALSE)
    }
    invisible(labels)
}
