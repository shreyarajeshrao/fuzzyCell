# The parent population each pair's `to` label is compared against in relative
# co-localisation: its own group, or its override. Pairs of a label with
# itself, and pairs whose `to` has no parent other than itself, are dropped.
# Returns the parent table spicyR::kontextualEngine() takes.
relativeTriples <- function(pairs, cellTypeTree, overrides) {
    cellMap <- buildCellTypeMap(cellTypeTree)$cellMap
    leafToGroup <- stats::setNames(cellMap$group, cellMap$cellType)
    parent <- ifelse(pairs$to %in% names(overrides),
        unname(overrides[pairs$to]), unname(leafToGroup[pairs$to]))
    keep <- !is.na(parent) & parent != pairs$to & pairs$from != pairs$to
    out <- data.frame(from = pairs$from[keep], to = pairs$to[keep],
        parent_name = parent[keep])
    parentLeaves <- labelLeaves(cellTypeTree, unique(out$parent_name))
    out$parent <- unname(parentLeaves[out$parent_name])
    out
}

# The cell types below each label (a cell type, group, wrapper or lineage), as
# a list named by label.
labelLeaves <- function(cellTypeTree, labels) {
    stats::setNames(lapply(labels, function(l) {
        resolveTreeLeaves(cellTypeTree, l)
    }), labels)
}

# Every ordered pair of labels the co-localisation could test, as three data
# frames of `from` and `to`: group x group, group x cell type (both
# directions) and cell type x cell type. The groups are those in confList.
generateComparisons <- function(cellTypeTree, confList) {
    groups <- names(confList)
    leaves <- buildCellTypeMap(cellTypeTree)$cellMap$cellType
    groupLeaf <- expand.grid(group = groups, leaf = leaves,
        stringsAsFactors = FALSE)
    list(
        groupGroup = expand.grid(from = groups, to = groups,
            stringsAsFactors = FALSE),
        groupLeaf = rbind(
            data.frame(from = groupLeaf$group, to = groupLeaf$leaf),
            data.frame(from = groupLeaf$leaf, to = groupLeaf$group)),
        leafLeaf = expand.grid(from = leaves, to = leaves,
            stringsAsFactors = FALSE)
    )
}

# The pairs of labels to test: every pair of groups and cell types, except a
# label with itself or with a group that contains it, as these are not two
# distinct cell populations.
colocalisationPairs <- function(cellTypeTree, confList) {
    comparisons <- generateComparisons(cellTypeTree, confList)
    pairs <- unique(do.call(rbind, lapply(comparisons, function(d) {
        data.frame(from = d$from, to = d$to)
    })))
    leaves <- labelLeaves(cellTypeTree, unique(c(pairs$from, pairs$to)))
    covers <- function(a, b) {
        vapply(seq_along(a), function(k) b[k] %in% leaves[[a[k]]], TRUE)
    }
    keep <- pairs$from != pairs$to & !covers(pairs$from, pairs$to) &
        !covers(pairs$to, pairs$from)
    out <- pairs[keep, , drop = FALSE]
    rownames(out) <- NULL
    out
}

# The rows of `tab` whose from and to form one of `pairs`, or all of them if
# `pairs` is NULL.
keepPairs <- function(tab, pairs) {
    if (is.null(pairs)) return(tab)
    wanted <- paste(pairs$from, pairs$to, sep = "__")
    tab[paste(tab$from, tab$to, sep = "__") %in% wanted, , drop = FALSE]
}

# Per-image effects of a spicyR result in long format (imageID, from, to,
# colocalisation, and parent for relative co-localisation).
perImageEffects <- function(res) {
    comp <- res$comparisons
    out <- do.call(rbind, lapply(seq_len(nrow(comp)), function(k) {
        d <- data.frame(imageID = res$imageIDs, from = comp$from[k],
            to = comp$to[k],
            colocalisation = res$pairwiseAssoc[[comp$labels[k]]])
        if (!is.null(comp$parent)) d$parent <- comp$parent[k]
        d
    }))
    out[!is.na(out$colocalisation), , drop = FALSE]
}

# The model table of a spicyR result, for the requested pairs: the estimate
# (second level of the outcome minus the first), its standard error,
# statistic, p-value and BH-adjusted p-value. Takes a spicy(method = "image")
# result or a kontextualEngine() result.
modelTable <- function(res, pairs = NULL) {
    if (is.null(res$cellResults)) return(imageModelTable(res, pairs))
    tab <- keepPairs(res$cellResults, pairs)
    out <- data.frame(from = tab$from, to = tab$to)
    if (!is.null(tab$parent)) out$parent <- tab$parent
    out$estimate <- tab$excess_difference
    out$se <- tab$se
    out$statistic <- tab$excess_difference / tab$se
    out$p.value <- tab$p_value
    if (!is.null(tab$side)) out$side <- tab$side
    out$p.adj <- stats::p.adjust(out$p.value, "BH")
    rownames(out) <- NULL
    out
}

# The model table of a spicy(method = "image") result: one linear model of the
# per-image co-localisation on the outcome per pair, whose coefficient for the
# second level of the outcome is the estimate.
imageModelTable <- function(res, pairs = NULL) {
    comp <- res$comparisons
    out <- data.frame(from = comp$from, to = comp$to,
        estimate = res$coefficient[[2]], se = res$se[[2]],
        statistic = res$statistic[[2]], p.value = res$p.value[[2]])
    out <- keepPairs(out, pairs)
    out$p.adj <- stats::p.adjust(out$p.value, "BH")
    rownames(out) <- NULL
    out
}

#' Co-localisation analysis with soft cell type weights
#'
#' Tests, for every pair of leaf cell types and groups, whether the
#' co-localisation of the two differs between the two conditions of an outcome,
#' once with the original labels and once with the bootstrap confidences as soft
#' cell type weights. The co-localisation of each image is the cross-type
#' L-function, computed by spicyR (`spicyR::spicy(method = "image")`), and
#' relative co-localisation is Statial's Kontextual test
#' (`spicyR::kontextualEngine()`).
#'
#' @details
#' Co-localisation is the cross-type L-function in each image: positive values
#' mean the `from` cells cluster around the `to` cells within radius
#' `r`, negative values that they avoid them. Each pair is tested with a
#' linear model of the per-image values on the outcome (a mixed model if some
#' subjects have several images).
#'
#' With weights, each cell counts towards each subtype of its group in
#' proportion to its bootstrap confidence (see
#' [confidenceWeights()]).
#'
#' Relative co-localisation compares the co-localisation of `from` with
#' `to` to that with a parent population of `to`: by default `to`'s
#' own group (see `parentOverrides`). Attraction that comes only from both
#' cell types preferring the same region is then not counted. Pairs whose `to`
#' has no parent are not tested, and no label is paired with itself or with a
#' group that contains it.
#'
#' @param spe A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a
#'   cell ID, cell type, image and coordinates for every cell (see
#'   [getCellData()]).
#' @param cellTypeTree A cell type tree, see [cellTypeTree].
#' @param confList Confidence list from
#'   [bootstrapLineageConfidence()].
#' @param outcomeCol Name of the `colData(spe)` column with the per-image
#'   outcome. It must have exactly two levels and a single value per image.
#' @param subject Name of the `colData(spe)` column identifying the subject
#'   (for example the patient) each image comes from. With `NULL` every image
#'   is its own subject.
#' @param r Radius within which cells count as neighbours, in the units of the
#'   coordinates; a single positive number.
#' @param window The area of each image that densities are computed over:
#'   `"convex"` (the default), the convex hull of the image's cells, or
#'   `"square"`, the rectangle around them. The same window is used for
#'   co-localisation and relative co-localisation, and
#'   [runNicheIdentification()] has the same
#'   default.
#' @param cores Number of cores spicyR runs on.
#' @param parentOverrides Named character vector mapping a `to` label to the
#'   ancestor node (a group, wrapper or lineage in `cellTypeTree`) that
#'   relative co-localisation should compare against instead of the label's own
#'   group, for example `c("Fibroblast lineage" = "Mesenchymal")`.
#' @param autoOverrideParents Whether to add overrides found by
#'   [autoParentOverrides()]. Use `FALSE` to
#'   condition only on the parents you list in `parentOverrides`.
#' @param cellType,cellID,imageID Names of the `colData(spe)` columns
#'   holding each cell's cell type label, a unique cell ID (matching the
#'   `cellID` column of `confList`), and the image it is in.
#' @param x,y Names of the columns holding the cell coordinates. For a
#'   `SpatialExperiment` whose coordinates are in `spatialCoords()`, the
#'   defaults `"x"` and `"y"` refer to its first two coordinates.
#'
#' @return A list with:
#' * the tests, one row per pair: `modelWeighted`, `modelUnweighted`,
#'   `modelRelativeWeighted` and `modelRelativeUnweighted`. Columns:
#'   `from`, `to`, `estimate` (the second level of the outcome minus
#'   the first), `se`, `statistic`, `p.value` and the BH-adjusted
#'   `p.adj`; the relative tables also have the `parent` and `side`
#'   (`"attract"` or `"avoid"`).
#' * the per-image values: `weighted`, `unweighted`,
#'   `relativeWeighted` and `relativeUnweighted`. Columns: `imageID`,
#'   `from`, `to`, `colocalisation`, and `parent` for the
#'   relative versions.
#' * `spicy`: the spicyR results themselves (`weighted`,
#'   `unweighted`, `relativeWeighted` and `relativeUnweighted`), for
#'   spicyR's plotting functions.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' # a small example: the endothelial cells and granulocytes only
#' groups <- c("Endothelial lineage", "Granulocytes")
#' tree <- lungCellTypeTree[["Immune"]][["Myeloid"]]["Granulocytes"]
#' tree <- c(lungCellTypeTree["Endothelial lineage"], tree)
#' cellTypes <- unlist(lapply(groups, function(g) resolveTreeLeaves(tree, g)))
#' spe <- lungSubset[, lungSubset$cellType %in% cellTypes]
#' res <- runColocalisationModels(spe, tree, lungConfList[groups],
#'     outcomeCol = "disease_status")
#' head(res$modelWeighted)
#' @export
runColocalisationModels <- function(spe, cellTypeTree, confList, outcomeCol,
        r = 50, window = c("convex", "square"), subject = NULL, cores = 1,
        parentOverrides = character(0),
        autoOverrideParents = TRUE, cellType = "cellType", cellID = "cellID",
        imageID = "imageID", x = "x", y = "y") {
    # check the arguments before doing any work
    window <- match.arg(window)
    checkSpe(spe)
    checkRadii(r)
    if (length(r) != 1) {
        stop("`r` must be a single positive number (the radius)", call. = FALSE)
    }
    if (!is.logical(autoOverrideParents) || length(autoOverrideParents) != 1 ||
        is.na(autoOverrideParents)) {
        stop("`autoOverrideParents` must be TRUE or FALSE", call. = FALSE)
    }

    cells <- standardiseCells(getCellData(spe), list(cellID = cellID,
        cellType = cellType, imageID = imageID, x = x, y = y))
    checkOutcome(cells, outcomeCol, subject)
    if (!is.null(subject)) cells$.subject <- cells[[subject]]

    message("Generating comparisons...")
    pairs <- colocalisationPairs(cellTypeTree, confList)
    overrides <- parentOverrides
    if (autoOverrideParents) {
        auto <- autoParentOverrides(cellTypeTree)
        overrides <- c(auto[!names(auto) %in% names(parentOverrides)],
            parentOverrides)
    }
    triples <- relativeTriples(pairs, cellTypeTree, overrides)

    spicyArgs <- list(cells = cells, outcomeCol = outcomeCol,
        subject = if (is.null(subject)) NULL else ".subject", r = r,
        window = window, pairs = pairs, triples = triples, cores = cores)
    weighted <- do.call(runSpicy, c(spicyArgs, list(label = "weighted",
        weights = weightMatrix(cells, cellTypeTree, confList, TRUE))))
    unweighted <- do.call(runSpicy, c(spicyArgs, list(label = "unweighted",
        weights = weightMatrix(cells, cellTypeTree, confList, FALSE))))
    res <- list(weighted = weighted$colocalisation,
        relativeWeighted = weighted$relative,
        unweighted = unweighted$colocalisation,
        relativeUnweighted = unweighted$relative)
    colocalisationResults(Filter(Negate(is.null), res), pairs)
}

# Runs spicyR's co-localisation test on `pairs`, and its relative
# co-localisation test on `triples`, with the cell weights `weights`. `label`
# ("weighted" or "unweighted") is used in the progress messages. Returns a list
# with `colocalisation` and `relative` (NULL when there are no triples).
runSpicy <- function(cells, outcomeCol, subject, r, window, pairs, triples,
        cores, weights, label) {
    message("Testing ", label, " co-localisation...")
    # the image method tests from[i] against to[i]; images are not weighted by
    # their number of cells
    colocalisation <- spicyR::spicy(cells, condition = outcomeCol,
        subject = subject, r = r, from = pairs$from, to = pairs$to,
        method = "image", window = window, weights = FALSE, cores = cores,
        cellWeights = weights)
    message("Testing ", label, " relative co-localisation...")
    # kontextualEngine() calls the rectangle around the cells "rectangle"
    relative <- if (nrow(triples) > 0) {
        spicyR::kontextualEngine(cells, triples, condition = outcomeCol,
            subject = subject, r = r, cellWeights = weights,
            window = if (window == "square") "rectangle" else "convex")
    }
    list(colocalisation = colocalisation, relative = relative)
}

# The result list of runColocalisationModels() from the spicyR results `res`.
colocalisationResults <- function(res, pairs) {
    # if no pair has a parent, the relative tables are empty but keep their
    # columns
    noModel <- data.frame(from = character(), to = character(),
        parent = character(), estimate = numeric(), se = numeric(),
        statistic = numeric(), p.value = numeric(), p.adj = numeric())
    noImages <- data.frame(imageID = character(), from = character(),
        to = character(), colocalisation = numeric(), parent = character())
    model <- function(x, p = NULL) if (is.null(x)) noModel else modelTable(x, p)
    perImage <- function(x) if (is.null(x)) noImages else perImageEffects(x)
    list(
        weighted = perImage(res$weighted),
        unweighted = perImage(res$unweighted),
        modelWeighted = model(res$weighted, pairs),
        modelUnweighted = model(res$unweighted, pairs),
        relativeWeighted = perImage(res$relativeWeighted),
        relativeUnweighted = perImage(res$relativeUnweighted),
        modelRelativeWeighted = model(res$relativeWeighted),
        modelRelativeUnweighted = model(res$relativeUnweighted),
        spicy = res
    )
}
