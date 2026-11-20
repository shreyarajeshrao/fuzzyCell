# The proportion of each image's cells in each value of the column `column` of
# `cells` (cell types or groups): an images x values table whose rows sum to 1.
imageProportions <- function(cells, column) {
    tab <- table(sample = cells$imageID, clusters = cells[[column]])
    tab / rowSums(tab)
}

# Every subtype of every group in confList (except skipGroups), as a data
# frame with one row per (group, subtype).
groupSubtypes <- function(confList, skipGroups = character(0)) {
    groups <- setdiff(names(confList), skipGroups)
    do.call(rbind, lapply(groups, function(group) {
        data.frame(group = group,
            subtype = setdiff(colnames(confList[[group]]), "cellID"))
    }))
}

# A model of each column of `props` (images x cell types) on the outcome, one
# per cell type: a linear model, or with `mixed = TRUE` a linear mixed model
# with a random intercept per subject (mixedModelTest()). `estimate` is the
# second level of the outcome minus the first. The outcome and subject are
# checked beforehand, in proportionCells().
runProportionModels <- function(props, cells, outcomeCol, subject = NULL,
        mixed = FALSE, label = "cell type proportions") {
    keep <- c("imageID", outcomeCol, if (mixed) subject)
    outcomeDat <- unique(cells[, keep, drop = FALSE])
    propDat <- data.frame(imageID = rownames(props), unclass(props),
        check.names = FALSE)
    rows <- match(propDat$imageID, outcomeDat$imageID)
    propDat[[outcomeCol]] <- outcomeDat[[outcomeCol]][rows]
    if (mixed) propDat$.subject <- outcomeDat[[subject]][rows]
    noModel <- function(cellType) {
        data.frame(cellType = cellType, estimate = NA, statistic = NA,
            p.value = NA, singular = NA)
    }

    out <- do.call(rbind, lapply(colnames(props), function(cellType) {
        dat <- propDat[!is.na(propDat[[cellType]]) &
            !is.na(propDat[[outcomeCol]]), , drop = FALSE]
        # skip cell types with too little data to fit a meaningful model
        if (length(unique(dat[[outcomeCol]])) < 2 || nrow(dat) < 5) {
            return(noModel(cellType))
        }
        if (mixed) {
            fit <- mixedModelTest(dat, cellType, outcomeCol)
            if (is.null(fit)) return(noModel(cellType))
            return(data.frame(cellType = cellType,
                estimate = fit[["estimate"]], statistic = fit[["statistic"]],
                p.value = fit[["p.value"]], singular = fit[["singular"]]))
        }
        formula <- as.formula(paste0("`", cellType, "` ~ `", outcomeCol, "`"))
        fit <- tryCatch(lm(formula, data = dat), error = function(e) NULL)
        if (is.null(fit)) return(noModel(cellType))
        coef <- summary(fit)$coefficients
        data.frame(cellType = cellType, estimate = coef[2, "Estimate"],
            statistic = coef[2, "t value"], p.value = coef[2, "Pr(>|t|)"],
            singular = NA)
    }))
    out <- reportSingular(out, label)
    out$p.adj <- p.adjust(out$p.value, method = "BH")
    rownames(out) <- NULL
    out
}

# Whether the models are mixed models: TRUE when `subject` is given and some
# subject has several images, which needs lmerTest; a message says so.
useMixedModels <- function(cells, subject) {
    if (!hasRepeatedImages(cells, subject)) return(FALSE)
    needLmerTest()
    message("Some subjects have several images: fitting linear mixed models ",
        "with a random intercept per subject.")
    TRUE
}

# A proportions table (images x cell types) in long format: imageID, cellType
# and the proportion in a column named `valueCol`.
propsToLong <- function(props, valueCol) {
    out <- data.frame(
        imageID = rep(rownames(props), times = ncol(props)),
        cellType = rep(colnames(props), each = nrow(props)),
        value = as.vector(props)
    )
    names(out)[3] <- valueCol
    out
}

# Weighted proportion of every subtype of every group in confList, per image:
# the sum of the confidences of the group's cells for that subtype, divided by
# all cells in the image.
runWeightedProportions <- function(cells, cellTypeTree, confList,
        skipGroups = character(0)) {
    cells <- annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)

    # these do not depend on the group or subtype, so they are computed once
    # instead of once per (group, subtype) pair
    trueTotal <- table(cells$imageID)
    allImages <- names(trueTotal)
    rowsByGroup <- split(seq_len(nrow(cells)), cells$group)
    pairs <- groupSubtypes(confList, skipGroups)

    out <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
        group <- pairs$group[i]
        subtype <- pairs$subtype[i]

        # every member of the group counts towards `subtype` with weight equal
        # to its bootstrap confidence of actually being that subtype
        inGroupRows <- rowsByGroup[[group]]
        if (is.null(inGroupRows)) inGroupRows <- integer(0)
        conf <- confList[[group]]
        groupWeight <- conf[[subtype]][match(cells$cellID[inGroupRows],
            conf$cellID)]
        subtypeSum <- tapply(groupWeight, cells$imageID[inGroupRows], sum)

        full <- setNames(rep(0, length(allImages)), allImages)
        full[names(subtypeSum)] <- subtypeSum
        data.frame(imageID = allImages, group = group, cellType = subtype,
            proportionWeighted = unname(full /
                as.vector(trueTotal[allImages])))
    }))
    rownames(out) <- NULL
    out
}

# Checks the arguments the proportion pipelines share: `spe` is a
# SummarizedExperiment, `confList` has the structure that
# bootstrapLineageConfidence() returns for groups of the tree, and `skipGroups`
# names groups of the tree (the tree itself, not confList, so the same
# skipGroups can be passed here and to bootstrapLineageConfidence()).
checkProportionArgs <- function(spe, cellTypeTree, confList, skipGroups) {
    checkSpe(spe)
    # also validates the tree
    groups <- buildCellTypeMap(cellTypeTree)$groupsWithChildren
    shownGroups <- paste(utils::head(groups, 10), collapse = ", ")
    if (!is.list(confList) || is.data.frame(confList) ||
        length(confList) == 0 || is.null(names(confList)) ||
        any(!nzchar(names(confList)))) {
        stop("`confList` must be a named list with one data.frame per group, ",
            "as returned by bootstrapLineageConfidence()", call. = FALSE)
    }
    # a loop, so that the first malformed group stops the check
    for (g in names(confList)) {
        conf <- confList[[g]]
        if (!is.data.frame(conf) || !"cellID" %in% names(conf) ||
            ncol(conf) < 2) {
            stop("`confList[[\"", g, "\"]]` must be a data.frame with a ",
                "`cellID` column and one column of confidence scores per ",
                "subtype, as returned by bootstrapLineageConfidence()",
                call. = FALSE)
        }
    }
    checkNamesAreGroups(names(confList), groups, "`confList` has",
        "group", shownGroups)
    if (!is.character(skipGroups) || anyNA(skipGroups)) {
        stop("`skipGroups` must be a character vector of group names (use ",
            "character(0) to skip none)", call. = FALSE)
    }
    checkNamesAreGroups(skipGroups, groups, "`skipGroups` has", "name",
        shownGroups)
    if (length(setdiff(names(confList), skipGroups)) == 0) {
        stop("Every group in `confList` is in `skipGroups`, so there is ",
            "nothing left to analyse", call. = FALSE)
    }
    invisible(TRUE)
}

# Stops if any of `labels` is not a group of the tree: "<what> 2 <noun>s that
# are not a group in cellTypeTree: ...".
checkNamesAreGroups <- function(labels, groups, what, noun, shownGroups) {
    unknown <- setdiff(labels, groups)
    if (length(unknown) > 0) {
        shown <- paste0("'", utils::head(unknown, 5), "'", collapse = ", ")
        stop(what, " ", length(unknown), " ", noun,
            if (length(unknown) > 1) "s", " that ",
            if (length(unknown) > 1) "are" else "is", " not a group in ",
            "cellTypeTree: ", shown, ". Groups in the tree: ", shownGroups,
            call. = FALSE)
    }
    invisible(TRUE)
}

# The per-cell data both proportion pipelines start from, after checking their
# arguments: the standard columns, each cell's group and lineage from the tree
# (an ungrouped cell type is its own group), and a checked outcome and subject.
proportionCells <- function(spe, cellTypeTree, confList, outcomeCol, skipGroups,
        subject, cellType, cellID, imageID) {
    checkProportionArgs(spe, cellTypeTree, confList, skipGroups)
    cells <- standardiseCells(getCellData(spe), list(cellID = cellID,
        cellType = cellType, imageID = imageID))
    cells <- annotateCellTypeTree(cells, cellTypeTree)
    cells$group[is.na(cells$group)] <- cells$cellType[is.na(cells$group)]
    checkOutcome(cells, outcomeCol, subject)
    cells
}

#' Cell type proportion analysis
#'
#' Computes per-image proportions of every cell type using original labels
#' (unweighted) and bootstrap confidences (weighted), and tests each against an
#' outcome. `runProportionsRelativeToParentPipeline()` instead expresses each
#' subtype as a proportion of its parent group.
#'
#' @details
#' In `runProportionsPipeline()` the unweighted results cover every cell type
#' and also the whole groups that were bootstrapped. Only the subtypes of those
#' groups have confidences, so for every other cell type (cell types without a
#' group, and the group-level proportions) the weighted proportion is the
#' unweighted one: their labels are not uncertain. The weighted and unweighted
#' models therefore test the same set of cell types, and their adjusted p-values
#' are comparable.
#'
#' Each cell type's proportion is tested with a linear model on the outcome,
#' with each image as one observation. When `subject` is given and some
#' subject has several images, those images are not independent, so a linear
#' mixed model with a random intercept per subject is fitted instead
#' (`lmerTest::lmer()`, as spicyR fits for co-localisation); this needs the
#' lmerTest package, and a message says so.
#'
#' @param spe A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a
#'   cell ID, cell type and image for every cell in `colData`.
#' @param cellTypeTree A cell type tree, see [cellTypeTree].
#' @param confList Confidence list from
#'   [bootstrapLineageConfidence()]: a
#'   named list with one data frame per group, whose names must be groups of
#'   `cellTypeTree`.
#' @param outcomeCol Name of the per-image outcome column in `colData(spe)`.
#'   It must have exactly two levels (subset your data to two conditions first)
#'   and a single value per image. An error is raised otherwise, before any
#'   computation.
#' @param skipGroups Names of groups of `cellTypeTree` to exclude, as a
#'   character vector (the default `character(0)` excludes none). An error is
#'   raised for a name that is not a group of the tree, so the same
#'   `skipGroups` can be used here and in
#'   [bootstrapLineageConfidence()].
#' @param subject Name of the `colData(spe)` column identifying the subject
#'   (for example the patient) each image comes from, or `NULL`. If given,
#'   every image must belong to one subject; subjects with several images are
#'   handled with mixed models (see Details). With `NULL` every image is
#'   treated as a separate subject, and a message says so.
#' @param cellType,cellID,imageID Names of the `colData(spe)` columns
#'   holding each cell's cell type label, a unique cell ID (matching the
#'   `cellID` column of `confList`), and the image it is in.
#'
#' @return
#' `runProportionsPipeline()` returns a list with:
#' * `unweighted`: the proportion of each image's cells in each cell type
#'   and bootstrapped group, from the original labels (columns `imageID`,
#'   `cellType`, `proportionUnweighted`);
#' * `weighted`: the same from the bootstrap confidences (columns
#'   `imageID`, `group`, `cellType`, `proportionWeighted`);
#' * `proportions`: both together (columns `imageID`,
#'   `cellType`, `proportionUnweighted`, `group`,
#'   `proportionWeighted`);
#' * `modelUnweighted`, `modelWeighted`: one row per cell type, with
#'   the `estimate` (the second level of the outcome minus the first), its
#'   `statistic` (t value), `p.value`, the BH-adjusted `p.adj` and the
#'   `group`.
#'
#' `runProportionsRelativeToParentPipeline()` returns `proportions` (each
#' subtype as a proportion of its group's cells, columns `imageID`,
#' `group`, `cellType`, `proportionUnweighted`,
#' `proportionWeighted`; `NaN` for images without cells of the group),
#' `modelUnweighted` and `modelWeighted`.
#'
#' In every table, `group` is the group a subtype belongs to; for a cell type
#' without a group, and for a whole group, it is the cell type or group itself.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' res <- runProportionsPipeline(lungSubset, lungCellTypeTree, lungConfList,
#'     outcomeCol = "disease_status", subject = "subject")
#' head(res$modelWeighted)
#' @export
runProportionsPipeline <- function(spe, cellTypeTree, confList, outcomeCol,
        skipGroups = character(0), subject = NULL, cellType = "cellType",
        cellID = "cellID", imageID = "imageID") {
    cells <- proportionCells(spe, cellTypeTree, confList, outcomeCol,
        skipGroups, subject, cellType, cellID, imageID)
    mixed <- useMixedModels(cells, subject)

    message("Computing unweighted proportions...")
    # the cell types, and the whole groups that were bootstrapped (an
    # ungrouped cell type is its own group, so its group column would only
    # repeat the cell type column)
    unweightedProps <- imageProportions(cells, "cellType")
    groupProps <- imageProportions(cells, "group")[rownames(unweightedProps),
        setdiff(names(confList), skipGroups), drop = FALSE]
    unweightedProps <- cbind(unweightedProps, groupProps)
    unweightedLong <- propsToLong(unweightedProps, "proportionUnweighted")
    modelUnweighted <- runProportionModels(unweightedProps, cells, outcomeCol,
        subject, mixed, "unweighted cell type proportions")

    message("Computing weighted proportions...")
    weighted <- withHardLabelled(runWeightedProportions(cells, cellTypeTree,
        confList, skipGroups = skipGroups), unweightedLong)

    message("Fitting weighted models...")
    weightedProps <- weighted |>
        tidyr::pivot_wider(id_cols = "imageID", names_from = "cellType",
            values_from = "proportionWeighted") |>
        tibble::column_to_rownames("imageID") |>
        as.matrix()
    modelWeighted <- runProportionModels(weightedProps, cells, outcomeCol,
        subject, mixed, "weighted cell type proportions")
    groupOf <- function(cellType) {
        weighted$group[match(cellType, weighted$cellType)]
    }
    modelUnweighted$group <- groupOf(modelUnweighted$cellType)
    modelWeighted$group <- groupOf(modelWeighted$cellType)

    list(
        proportions = joinProportions(unweightedLong, weighted),
        weighted = weighted,
        unweighted = unweightedLong,
        modelWeighted = modelWeighted,
        modelUnweighted = modelUnweighted
    )
}

# The weighted proportions `weighted` with every other cell type of
# `unweightedLong` added. Only subtypes of bootstrapped groups have
# confidences; everything else (ungrouped cell types and whole groups) has hard
# labels, so its weighted proportion is the unweighted one. This puts both sets
# of models on the same cell types, so their adjusted p-values are comparable.
withHardLabelled <- function(weighted, unweightedLong) {
    hard <- unweightedLong[!unweightedLong$cellType %in% weighted$cellType, ]
    out <- rbind(weighted, data.frame(imageID = hard$imageID,
        group = hard$cellType, cellType = hard$cellType,
        proportionWeighted = hard$proportionUnweighted))
    rownames(out) <- NULL
    out
}

# The unweighted and weighted proportions side by side, one row per image and
# cell type found in both.
joinProportions <- function(unweightedLong, weighted) {
    rows <- match(paste(unweightedLong$imageID, unweightedLong$cellType,
        sep = "__"), paste(weighted$imageID, weighted$cellType, sep = "__"))
    found <- !is.na(rows)
    out <- cbind(unweightedLong[found, ],
        weighted[rows[found], c("group", "proportionWeighted")])
    rownames(out) <- NULL
    out
}

# Proportion of each subtype within its parent group, per image, both
# unweighted (original labels) and weighted (bootstrap confidences). A long
# data.frame with the subtype's group in `group`; images with no cells of the
# group get NaN.
calcProportionsRelativeToParent <- function(cells, cellTypeTree, confList,
        skipGroups = character(0)) {
    cells <- annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    allImages <- unique(cells$imageID)
    rowsByGroup <- split(seq_len(nrow(cells)), cells$group)
    pairs <- groupSubtypes(confList, skipGroups)

    out <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
        group <- pairs$group[i]
        subtype <- pairs$subtype[i]
        inGroupRows <- rowsByGroup[[group]]
        if (is.null(inGroupRows)) inGroupRows <- integer(0)
        groupImages <- cells$imageID[inGroupRows]
        perImage <- function(counts) {
            full <- setNames(rep(0, length(allImages)), allImages)
            full[names(counts)] <- counts
            unname(full)
        }

        # denominator: hard count of cells in the parent group, per image
        parentTotal <- perImage(table(groupImages))
        # unweighted numerator: hard count of cells labelled `subtype`
        unweightedCount <- perImage(table(groupImages[
            cells$cellType[inGroupRows] == subtype]))
        # weighted numerator: the summed confidences for `subtype`
        conf <- confList[[group]]
        weight <- conf[[subtype]][match(cells$cellID[inGroupRows], conf$cellID)]
        weightedSum <- perImage(tapply(weight, groupImages, sum))

        data.frame(imageID = allImages, group = group, cellType = subtype,
            proportionUnweighted = unweightedCount / parentTotal,
            proportionWeighted = weightedSum / parentTotal)
    }))
    rownames(out) <- NULL
    out
}

#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' rel <- runProportionsRelativeToParentPipeline(lungSubset, lungCellTypeTree,
#'     lungConfList, outcomeCol = "disease_status", subject = "subject")
#' head(rel$modelWeighted)
#' @rdname runProportionsPipeline
#' @export
runProportionsRelativeToParentPipeline <- function(spe, cellTypeTree, confList,
        outcomeCol, skipGroups = character(0), subject = NULL,
        cellType = "cellType", cellID = "cellID", imageID = "imageID") {
    cells <- proportionCells(spe, cellTypeTree, confList, outcomeCol,
        skipGroups, subject, cellType, cellID, imageID)
    mixed <- useMixedModels(cells, subject)

    message("Computing proportions relative to parent...")
    proportions <- calcProportionsRelativeToParent(cells, cellTypeTree,
        confList, skipGroups = skipGroups)
    groupOf <- function(cellType) {
        proportions$group[match(cellType, proportions$cellType)]
    }

    message("Fitting unweighted models...")
    unweightedProps <- proportions |>
        tidyr::pivot_wider(id_cols = "imageID", names_from = "cellType",
            values_from = "proportionUnweighted") |>
        tibble::column_to_rownames("imageID") |>
        as.matrix()
    modelUnweighted <- runProportionModels(unweightedProps, cells, outcomeCol,
        subject, mixed, "unweighted proportions relative to the parent")
    modelUnweighted$group <- groupOf(modelUnweighted$cellType)

    message("Fitting weighted models...")
    weightedProps <- proportions |>
        tidyr::pivot_wider(id_cols = "imageID", names_from = "cellType",
            values_from = "proportionWeighted") |>
        tibble::column_to_rownames("imageID") |>
        as.matrix()
    modelWeighted <- runProportionModels(weightedProps, cells, outcomeCol,
        subject, mixed, "weighted proportions relative to the parent")
    modelWeighted$group <- groupOf(modelWeighted$cellType)

    list(
        proportions = proportions,
        modelWeighted = modelWeighted,
        modelUnweighted = modelUnweighted
    )
}
