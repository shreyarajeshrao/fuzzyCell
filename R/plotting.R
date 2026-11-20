# Shared helpers ----

# Stops with an install hint when a Suggests package needed by `caller` is
# missing.
checkInteractivePackages <- function(caller) {
    packages <- c("ggiraph", "gdtools")
    missing <- packages[!vapply(packages, requireNamespace, logical(1),
        quietly = TRUE)]
    if (length(missing) > 0) {
        needed <- paste(missing, collapse = " and ")
        install <- toString(dQuote(missing, FALSE))
        stop(caller, "() draws an interactive plot, which needs the ", needed,
            " package", if (length(missing) > 1) "s", ". Install with ",
            "BiocManager::install(c(", install, ")).", call. = FALSE)
    }
    invisible(TRUE)
}

# Makes a ggplot interactive with ggiraph. By default girafe() embeds the font
# files it measured the text with, several MB in every HTML page; here the text
# is measured with Liberation Sans, which gdtools ships, and drawn by the
# browser with Liberation Sans or a metric-compatible system font, so nothing
# is embedded.
interactivePlot <- function(p) {
    fonts <- gdtools::font_set(sans = gdtools::font_liberation("sans"))
    fonts$dependencies <- list()
    widget <- ggiraph::girafe(ggobj = p, height_svg = 3.5, font_set = fonts)
    widget$x$html <- gsub("font-family='Liberation Sans'",
        "font-family='Liberation Sans, Arial, Helvetica, sans-serif'",
        widget$x$html, fixed = TRUE)
    widget
}

# The per-cell data of `spe`, after the checks every spatial plot shares:
# `cellSize`, `spe`, and that `imageID` is an image of the data.
spatialCells <- function(spe, imageID, cellSize) {
    checkCellSize(cellSize)
    checkSpe(spe)
    cells <- getCellData(spe)
    checkImage(cells, imageID)
    cells
}

# The cells of image `imageID`, with each cell's group and lineage from the
# tree if one is given.
imageSubset <- function(cells, imageID, cellTypeTree = NULL) {
    if (!is.null(cellTypeTree)) {
        cells <- annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    }
    cells[cells$imageID == imageID, , drop = FALSE]
}

# A spatial plot with the `background` cells drawn in grey; `...` (such as
# alpha) is passed to their layer.
spatialPlot <- function(background, ...) {
    ggplot() +
        geom_point(data = background, aes(x = .data$x, y = .data$y),
            colour = "grey80", size = 0.7, ...)
}

# The layers that finish every spatial plot: image coordinates (y pointing
# down), equal axis scales, the axis titles and the `title`; `...` are further
# labels.
spatialFinish <- function(title, ...) {
    list(scale_y_reverse(), coord_equal(),
        labs(title = title, x = "x", y = "y", ...), theme_classic())
}

# Argument checks for the plots ----

# Checks that `cellSize`, the size of the highlighted cells in a spatial plot,
# is a single positive number.
checkCellSize <- function(cellSize) {
    if (!is.numeric(cellSize) || length(cellSize) != 1 || is.na(cellSize) ||
        !is.finite(cellSize) || cellSize <= 0) {
        stop("`cellSize` must be a single positive number", call. = FALSE)
    }
    invisible(cellSize)
}

# Checks that `imageID`, the image to plot, is a single image of the data.
checkImage <- function(cells, imageID) {
    images <- unique(as.character(cells$imageID))
    if (!is.character(imageID) || length(imageID) != 1 || is.na(imageID)) {
        stop("`imageID` must be a single string", call. = FALSE)
    }
    if (!imageID %in% images) {
        shown <- paste(utils::head(images, 10), collapse = ", ")
        stop("`imageID = \"", imageID, "\"` is not an image in the data. ",
            "Images: ", shown, if (length(images) > 10) ", ...",
            call. = FALSE)
    }
    invisible(imageID)
}

# Checks that `gene`, the argument called `arg`, is a single string naming a
# row of the count matrix.
checkGene <- function(counts, gene, arg, assayName) {
    if (!is.character(gene) || length(gene) != 1 || is.na(gene)) {
        stop("`", arg, "` must be a single string naming a gene",
            call. = FALSE)
    }
    if (!gene %in% rownames(counts)) {
        stop("`", arg, " = \"", gene, "\"` is not a gene (a row) of the '",
            assayName, "' assay", call. = FALSE)
    }
    invisible(gene)
}

# The groups of the cell type tree that have cells in the data. Without a tree
# the `group` column of the data is used.
groupsInData <- function(cells, cellTypeTree = NULL) {
    if (is.null(cellTypeTree)) {
        return(unique(as.character(cells$group[!is.na(cells$group)])))
    }
    cellMap <- buildCellTypeMap(cellTypeTree)$cellMap
    unique(cellMap$group[cellMap$cellType %in% cells$cellType &
        !is.na(cellMap$group)])
}

# Checks that every group the user passed exists in the data. `groups` is a
# named list: the names are the arguments, the values what the user passed.
checkGroupsInData <- function(cells, cellTypeTree, groups) {
    available <- groupsInData(cells, cellTypeTree)
    shown <- paste(utils::head(available, 10), collapse = ", ")
    # a loop, so that the first unknown group stops the check
    for (arg in names(groups)) {
        g <- groups[[arg]]
        if (!is.character(g) || length(g) != 1 || is.na(g)) {
            stop("`", arg, "` must be a single group name", call. = FALSE)
        }
        if (!g %in% available) {
            stop("`", arg, " = \"", g, "\"` is not a group in the data. ",
                "Groups: ", shown, if (length(available) > 10) ", ...",
                call. = FALSE)
        }
    }
    invisible(TRUE)
}

# Checks that every label the user passed is a cell type or a group of the
# data.
checkLabelsInData <- function(cells, cellTypeTree, labels) {
    available <- union(unique(as.character(cells$cellType)),
        groupsInData(cells, cellTypeTree))
    # a loop, so that the first unknown label stops the check
    for (arg in names(labels)) {
        l <- labels[[arg]]
        if (!is.character(l) || length(l) != 1 || is.na(l)) {
            stop("`", arg, "` must be a single cell type or group name",
                call. = FALSE)
        }
        if (!l %in% available) {
            stop("`", arg, " = \"", l, "\"` is neither a cell type nor a ",
                "group in the data", call. = FALSE)
        }
    }
    invisible(TRUE)
}

# Checks that `subtype` (the argument called `arg`) has confidence scores for
# `group` in confList; the group itself stands for the whole group.
checkSubtype <- function(confList, group, subtype, arg) {
    if (!is.character(subtype) || length(subtype) != 1 || is.na(subtype)) {
        stop("`", arg, "` must be a single cell type name", call. = FALSE)
    }
    if (identical(group, subtype)) {
        return(invisible(TRUE))
    }
    conf <- confList[[group]]
    if (!is.data.frame(conf)) {
        stop("`confList` has no confidence scores for the group '", group,
            "'", call. = FALSE)
    }
    available <- setdiff(names(conf), "cellID")
    if (!subtype %in% available) {
        shown <- paste(available, collapse = ", ")
        stop("`", arg, " = \"", subtype, "\"` is not a subtype of '", group,
            "' in `confList`. Subtypes: ", shown, call. = FALSE)
    }
    invisible(TRUE)
}

# Niche labels for the cells: `niche` (the argument called `arg`) is either a
# column of the cell data or a vector of niches named by cell ID, such as
# `nicheWeighted` from runNicheIdentification(); a vector is added to the cells
# as the column ".niche". Returns the cells and the name of the niche column.
withNicheColumn <- function(cells, niche, arg) {
    if (is.null(names(niche))) {
        checkColumn(cells, niche, arg)
        return(list(cells = cells, col = niche))
    }
    if (!any(cells$cellID %in% names(niche))) {
        stop("`", arg, "` must be a column of the data or a vector of niches ",
            "named by cell ID, with at least some of the cell IDs in `spe`",
            call. = FALSE)
    }
    cells$.niche <- unname(niche)[match(cells$cellID, names(niche))]
    list(cells = cells, col = ".niche")
}

# Checks that `niche` is one of the niches in the column `nicheCol`, which the
# error message calls `label`.
checkNiche <- function(cells, niche, nicheCol, label = nicheCol) {
    values <- cells[[nicheCol]]
    niches <- sort(unique(as.character(values[!is.na(values)])))
    if (!is.character(niche) || length(niche) != 1 || is.na(niche)) {
        stop("`niche` must be a single niche name", call. = FALSE)
    }
    if (!niche %in% niches) {
        shown <- paste(niches, collapse = ", ")
        stop("`niche = \"", niche, "\"` is not a niche in `", label, "`. ",
            "Niches: ", shown, call. = FALSE)
    }
    invisible(niche)
}

# Checks the two groups and cell types of plotWeightedSpatialPair(): the
# groups differ and have cells in the data, and each cell type has confidences
# in its group.
checkPairGroups <- function(cells, cellTypeTree, confList, groupA, subtypeA,
        groupB, subtypeB) {
    if (identical(groupA, groupB)) {
        stop("`groupA` and `groupB` are both '", groupA, "'. This plot ",
            "compares cell types from two different groups; to show the ",
            "confidences within one group, use plotWeightedSpatial().",
            call. = FALSE)
    }
    checkGroupsInData(cells, cellTypeTree, list(groupA = groupA,
        groupB = groupB))
    checkSubtype(confList, groupA, subtypeA, "subtypeA")
    checkSubtype(confList, groupB, subtypeB, "subtypeB")
    invisible(TRUE)
}

# Cell type labels, confidence scores and markers ----

# Relabels every cell of `group` as `newLabel` in the `misclassType` column.
relabelGroup <- function(cells, group, newLabel) {
    if (is.null(cells$misclassType)) cells$misclassType <- cells$cellType
    cells$misclassType[cells$group == group] <- newLabel
    cells
}

# Each cell's weight for `subtype`: its confidence for the cells of `group`,
# and 1 for every other cell. When `subtype` is the group itself (the whole
# group, not one of its subtypes) every cell keeps weight 1, as group
# membership has no confidence.
getSubtypeWeight <- function(cells, confList, group, subtype) {
    weight <- rep(1, nrow(cells))
    if (identical(group, subtype)) return(weight)
    inGroup <- which(cells$group == group)
    conf <- confList[[group]]
    weight[inGroup] <- conf[[subtype]][match(cells$cellID[inGroup],
        conf$cellID)]
    weight
}

# Each cell's weight for the two cell types of a pair plot: groupA cells get
# their confidence for subtypeA, groupB cells their confidence for subtypeB,
# and every other cell 1.
twoGroupWeights <- function(cells, confList, groupA, subtypeA, groupB,
        subtypeB) {
    weightA <- getSubtypeWeight(cells, confList, groupA, subtypeA)
    weightB <- getSubtypeWeight(cells, confList, groupB, subtypeB)
    ifelse(cells$group == groupA, weightA,
        ifelse(cells$group == groupB, weightB, 1))
}

# A marker's values rescaled to 0-1 by their range within the cells `mask`, so
# that the transparency scale uses its full range whatever the marker's level;
# 0 everywhere if the marker does not vary there.
rescaleWithin <- function(x, mask) {
    rng <- range(x[mask], na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    pmin(pmax((x - rng[1]) / diff(rng), 0), 1)
}

#' Spatial plots of cell types and confidence scores
#'
#' Plots the cells of one image, highlighting a cell type or group, or showing
#' bootstrap confidence scores. Each plot shows the cell positions (`x`,
#' `y`) with the other cells of the image in grey.
#'
#' * `plotWeightedSpatialPair()` shows two cell types from two different
#'   groups, with every cell of the two groups shaded by its confidence of being
#'   that cell type.
#' * `plotHardSpatialPair()` shows two cell types (or groups) with their
#'   original labels.
#' * `plotWeightedSpatial()` shows the cells of one group, shaded by their
#'   confidence of being one of its cell types.
#' * `plotOriginalLabelSpatial()` shows the cells of one group, coloured
#'   by their original cell type labels.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   a cell ID, cell type, image and coordinates for every cell in the form
#'   given by [getCellData()] (`cellID`,
#'   `cellType`, `imageID`, `x` and `y`).
#' @param cellTypeTree A cell type tree, see [cellTypeTree]. For
#'   `plotHardSpatialPair()` it can be `NULL`: the tree is only needed to
#'   use a group name as a label.
#' @param confList Confidence list from
#'   [bootstrapLineageConfidence()].
#' @param imageID The image to plot. An error is raised if it is not an image
#'   in the data.
#' @param groupA,groupB,targetGroup Group names. An error is raised for a group
#'   that has no cells in the data.
#' @param subtypeA,subtypeB,targetSubtype Cell types within the groups, which
#'   must have confidence scores in `confList` (for
#'   `plotHardSpatialPair()`, a cell type or a group name in the data).
#' @param encoding How to show the confidence scores: `"alpha"` (the
#'   default), as transparency, or `"colour"`, as a colour gradient. Either way
#'   each of the two cell types has its own legend.
#' @param cellSize Size of the highlighted cells, the ones drawn in colour. The
#'   grey background cells always have size 0.7.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' plotHardSpatialPair(lungSubset, lungCellTypeTree, "VUHD113", "Mast",
#'     "Proliferating FBs")
#' plotWeightedSpatial(lungSubset, lungCellTypeTree, lungConfList, "VUHD113",
#'     "Endothelial lineage", "Arteriole")
#' plotOriginalLabelSpatial(lungSubset, lungCellTypeTree, "VUHD113",
#'     "Endothelial lineage")
#' plotWeightedSpatialPair(lungSubset, lungCellTypeTree, lungConfList,
#'     "VUHD113", "Granulocytes", "Mast", "Fibroblast lineage",
#'     "Proliferating FBs")
#' @name plotSpatial
#' @export
plotWeightedSpatialPair <- function(spe, cellTypeTree, confList, imageID,
        groupA, subtypeA, groupB, subtypeB, encoding = c("alpha", "colour"),
        cellSize = 1) {
    cells <- spatialCells(spe, imageID, cellSize)
    encoding <- match.arg(encoding)
    checkPairGroups(cells, cellTypeTree, confList, groupA, subtypeA, groupB,
        subtypeB)
    cells <- imageSubset(cells, imageID, cellTypeTree)
    cells <- relabelGroup(relabelGroup(cells, groupA, subtypeA), groupB,
        subtypeB)
    cells$weight <- twoGroupWeights(cells, confList, groupA, subtypeA, groupB,
        subtypeB)
    cellsA <- cells[cells$misclassType %in% subtypeA, , drop = FALSE]
    cellsB <- cells[cells$misclassType %in% subtypeB, , drop = FALSE]
    nameA <- paste0(subtypeA, "\nconfidence")
    nameB <- paste0(subtypeB, "\nconfidence")

    p <- spatialPlot(cells[!cells$misclassType %in% c(subtypeA, subtypeB), ,
        drop = FALSE], alpha = 0.5)
    # each cell type has its own scale and legend (ggnewscale)
    if (encoding == "alpha") {
        p <- p +
            geom_point(data = cellsA, aes(x = .data$x, y = .data$y,
                alpha = .data$weight), colour = "red", size = cellSize) +
            scale_alpha_continuous(limits = c(0, 1), name = nameA) +
            ggnewscale::new_scale("alpha") +
            geom_point(data = cellsB, aes(x = .data$x, y = .data$y,
                alpha = .data$weight), colour = "blue", size = cellSize) +
            scale_alpha_continuous(limits = c(0, 1), name = nameB)
    } else {
        p <- p +
            geom_point(data = cellsA, aes(x = .data$x, y = .data$y,
                colour = .data$weight), size = cellSize) +
            scale_colour_gradient(low = "grey80", high = "red",
                limits = c(0, 1), name = nameA) +
            ggnewscale::new_scale_colour() +
            geom_point(data = cellsB, aes(x = .data$x, y = .data$y,
                colour = .data$weight), size = cellSize) +
            scale_colour_gradient(low = "grey80", high = "blue",
                limits = c(0, 1), name = nameB)
    }
    p + spatialFinish(paste0(imageID, "\n", subtypeA, " (red) vs ", subtypeB,
        " (blue)"))
}

#' @rdname plotSpatial
#' @export
plotHardSpatialPair <- function(spe, cellTypeTree, imageID, subtypeA,
        subtypeB, cellSize = 1) {
    cells <- spatialCells(spe, imageID, cellSize)
    checkLabelsInData(cells, cellTypeTree, list(subtypeA = subtypeA,
        subtypeB = subtypeB))
    # with a tree the cells get a `group` column, so that a group name can be
    # used as a label; without one, a group name matches nothing
    cells <- imageSubset(cells, imageID, cellTypeTree)

    # a label is a whole group if it is one, and otherwise a cell type
    matchLabel <- function(label) {
        if (!is.null(cells$group) && label %in% cells$group) {
            cells$group %in% label
        } else {
            cells$cellType %in% label
        }
    }
    matchA <- matchLabel(subtypeA)
    matchB <- matchLabel(subtypeB)

    spatialPlot(cells[!matchA & !matchB, , drop = FALSE], alpha = 0.5) +
        geom_point(data = cells[matchA, , drop = FALSE],
            aes(x = .data$x, y = .data$y), colour = "red", size = cellSize) +
        geom_point(data = cells[matchB, , drop = FALSE],
            aes(x = .data$x, y = .data$y), colour = "blue", size = cellSize) +
        spatialFinish(paste0(imageID, "\n", subtypeA, " (red) vs ", subtypeB,
            " (blue)"))
}

#' @rdname plotSpatial
#' @export
plotWeightedSpatial <- function(spe, cellTypeTree, confList, imageID,
        targetGroup, targetSubtype, cellSize = 1) {
    cells <- spatialCells(spe, imageID, cellSize)
    checkGroupsInData(cells, cellTypeTree, list(targetGroup = targetGroup))
    checkSubtype(confList, targetGroup, targetSubtype, "targetSubtype")
    cells <- imageSubset(cells, imageID, cellTypeTree)
    cells$weight <- getSubtypeWeight(cells, confList, targetGroup,
        targetSubtype)
    inGroup <- cells$group %in% targetGroup

    spatialPlot(cells[!inGroup, , drop = FALSE]) +
        geom_point(data = cells[inGroup, , drop = FALSE],
            aes(x = .data$x, y = .data$y, alpha = .data$weight),
            colour = "red", size = cellSize) +
        scale_alpha_continuous(range = c(0.05, 1), limits = c(0, 1)) +
        spatialFinish(paste0(imageID, " - ", targetSubtype, " confidence"),
            alpha = "Confidence")
}

#' @rdname plotSpatial
#' @export
plotOriginalLabelSpatial <- function(spe, cellTypeTree, imageID, targetGroup,
        cellSize = 1) {
    cells <- spatialCells(spe, imageID, cellSize)
    checkGroupsInData(cells, cellTypeTree, list(targetGroup = targetGroup))
    cells <- imageSubset(cells, imageID, cellTypeTree)
    inGroup <- cells$group %in% targetGroup

    spatialPlot(cells[!inGroup, , drop = FALSE]) +
        geom_point(data = cells[inGroup, , drop = FALSE],
            aes(x = .data$x, y = .data$y, colour = .data$cellType),
            size = cellSize) +
        spatialFinish(paste0(imageID, " - ", targetGroup,
            " (original labels)"), colour = "Cell type")
}

#' Spatial plots of marker expression within a group
#'
#' Plots one image, with the cells of a group shaded by the expression of one
#' marker (`plotOneMarkerSpatial()`), or the cells of two groups each shaded
#' by the expression of its own marker (`plotTwoMarkerSpatial()`). Expression
#' is normalised by each cell's library size, then rescaled to 0-1 within the
#' group, so that the shading uses its full range.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   a cell ID, cell type, image and coordinates for every cell in the form
#'   given by [getCellData()], and a counts assay.
#' @param cellTypeTree A cell type tree, used to find the cells of each group,
#'   or `NULL` if `colData(spe)` already has a `group` column.
#' @param imageID The image to plot.
#' @param group,groupA,groupB Group name(s) whose cells are shaded.
#' @param marker,markerA,markerB Gene name(s), which must be rows of the assay.
#' @param colour,colourA,colourB Colour(s) used for the shading.
#' @param cellSize Size of the highlighted cells, the ones drawn in colour. The
#'   grey background cells always have size 0.7.
#' @param assayName Name of the assay with the counts.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' plotOneMarkerSpatial(lungSubset, lungCellTypeTree, "VUHD113",
#'     "Endothelial lineage", "HEY1")
#' plotTwoMarkerSpatial(lungSubset, lungCellTypeTree, "VUHD113",
#'     "Endothelial lineage", "PECAM1", "purple", "Fibroblast lineage",
#'     "COL1A1", "green")
#' @name plotMarkerSpatial
#' @export
plotOneMarkerSpatial <- function(spe, cellTypeTree, imageID, group,
        marker, colour = "purple", cellSize = 1, assayName = "counts") {
    cells <- markerCells(spe, cellTypeTree, imageID, list(group = group),
        list(marker = marker), cellSize, assayName)
    inGroup <- cells$group %in% group
    cells$normScaled <- rescaleWithin(cells$marker, inGroup)

    spatialPlot(cells[!inGroup, , drop = FALSE]) +
        geom_point(data = cells[inGroup, , drop = FALSE],
            aes(x = .data$x, y = .data$y, alpha = .data$normScaled),
            colour = colour, size = cellSize) +
        scale_alpha_continuous(range = c(0.05, 1),
            name = paste0(marker, " (norm.)")) +
        spatialFinish(paste0(imageID, " - ", marker, " in ", group))
}

#' @rdname plotMarkerSpatial
#' @export
plotTwoMarkerSpatial <- function(spe, cellTypeTree, imageID, groupA,
        markerA, colourA = "purple", groupB, markerB, colourB = "green",
        cellSize = 1, assayName = "counts") {
    cells <- markerCells(spe, cellTypeTree, imageID,
        list(groupA = groupA, groupB = groupB),
        list(markerA = markerA, markerB = markerB), cellSize, assayName)
    inA <- cells$group %in% groupA
    inB <- cells$group %in% groupB
    cells$normAScaled <- rescaleWithin(cells$markerA, inA)
    cells$normBScaled <- rescaleWithin(cells$markerB, inB)

    spatialPlot(cells[!inA & !inB, , drop = FALSE]) +
        geom_point(data = cells[inA, , drop = FALSE],
            aes(x = .data$x, y = .data$y, alpha = .data$normAScaled),
            colour = colourA, size = cellSize) +
        geom_point(data = cells[inB, , drop = FALSE],
            aes(x = .data$x, y = .data$y, alpha = .data$normBScaled),
            colour = colourB, size = cellSize) +
        scale_alpha_continuous(range = c(0.05, 1),
            name = "Expression (norm.)") +
        spatialFinish(paste0(imageID, "\n", markerA, " in ", groupA, " (",
            colourA, ") / ", markerB, " in ", groupB, " (", colourB, ")"))
}

# The cells of image `imageID` for the marker plots, after checking the
# arguments: each cell's group, and for every marker in `markers` (a named list,
# its names being both the argument names and the new columns) the marker's
# counts divided by the cell's library size.
markerCells <- function(spe, cellTypeTree, imageID, groups, markers, cellSize,
        assayName) {
    cells <- spatialCells(spe, imageID, cellSize)
    checkAssay(spe, assayName)
    if (!is.null(cellTypeTree)) {
        cells <- annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    }
    if (is.null(cells$group)) {
        stop("`cellTypeTree` is needed to find the cells of the group",
            if (length(groups) > 1) "s", ", because colData(spe) has no ",
            "`group` column", call. = FALSE)
    }
    checkGroupsInData(cells, cellTypeTree, groups)
    countMat <- SummarizedExperiment::assay(spe, assayName)
    librarySize <- Matrix::colSums(countMat) # dense and sparse counts
    cells[names(markers)] <- lapply(names(markers), function(arg) {
        checkGene(countMat, markers[[arg]], arg, assayName)
        countMat[markers[[arg]], ] / librarySize
    })
    cells[cells$imageID == imageID, , drop = FALSE]
}

# Proportions ----
# Plot of the results of runProportionsPipeline() and
# runProportionsRelativeToParentPipeline().

#' Plot one cell type's proportion by outcome
#'
#' Shows, for one cell type (or group), each image's proportion measured with
#' the original labels (unweighted) and with the bootstrap confidence weights
#' (weighted), for each level of the outcome. For the results of
#' [runProportionsPipeline()] the proportion
#' is of all cells in the image; for the results of
#' [runProportionsRelativeToParentPipeline()]
#' it is of the cells of the parent group. Each point is one image.
#'
#' @param results Output of
#'   [runProportionsPipeline()] or
#'   [runProportionsRelativeToParentPipeline()].
#' @param cellType The cell type to show: one of the `cellType` values in
#'   `results$proportions`. For
#'   [runProportionsPipeline()] results this
#'   can also be one of the bootstrapped groups.
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   `imageID` and the outcome in `colData`.
#' @param outcomeCol Name of the outcome column in `colData(spe)`.
#'
#' @return An interactive plot, a `girafe` widget that needs the ggiraph
#'   package.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' res <- runProportionsPipeline(lungSubset, lungCellTypeTree, lungConfList,
#'     outcomeCol = "disease_status", subject = "subject")
#' rel <- runProportionsRelativeToParentPipeline(lungSubset, lungCellTypeTree,
#'     lungConfList, outcomeCol = "disease_status", subject = "subject")
#' if (requireNamespace("ggiraph", quietly = TRUE) &&
#'     requireNamespace("gdtools", quietly = TRUE)) {
#'     plotProportionBoxplot(res, "Arteriole", lungSubset, "disease_status")
#'     plotProportionBoxplot(rel, "Arteriole", lungSubset, "disease_status")
#' }
#' @export
plotProportionBoxplot <- function(results, cellType, spe, outcomeCol) {
    checkInteractivePackages("plotProportionBoxplot")
    checkSpe(spe)
    p <- proportionBoxplotGG(results, cellType, getCellData(spe), outcomeCol)
    interactivePlot(p)
}

# The per-image unweighted and weighted values of `dat` (columns imageID,
# `unweightedCol` and `weightedCol`) in long format, two rows per image (type
# "Unweighted" and "Weighted", value in plotValue), with each image's outcome.
weightedLong <- function(dat, unweightedCol, weightedCol, cells, outcomeCol) {
    outcomeDat <- unique(cells[, c("imageID", outcomeCol)])
    long <- data.frame(
        imageID = rep(as.character(dat$imageID), each = 2),
        type = rep(c("Unweighted", "Weighted"), times = nrow(dat)),
        plotValue = as.vector(rbind(dat[[unweightedCol]], dat[[weightedCol]]))
    )
    long[[outcomeCol]] <- outcomeDat[[outcomeCol]][match(long$imageID,
        outcomeDat$imageID)]
    long
}

# the ggplot behind plotProportionBoxplot(), before it is made interactive
proportionBoxplotGG <- function(results, cellType, cells, outcomeCol) {
    needed <- c("imageID", "cellType", "group", "proportionUnweighted",
        "proportionWeighted")
    if (!is.list(results) || !is.data.frame(results$proportions) ||
        !all(needed %in% names(results$proportions))) {
        stop("`results` must be the output of runProportionsPipeline() or ",
            "runProportionsRelativeToParentPipeline()", call. = FALSE)
    }
    available <- unique(results$proportions$cellType)
    if (!is.character(cellType) || length(cellType) != 1 || is.na(cellType)) {
        stop("`cellType` must be a single string", call. = FALSE)
    }
    if (!cellType %in% available) {
        shown <- paste(utils::head(available, 10), collapse = ", ")
        stop("`cellType = \"", cellType, "\"` is not in the results. ",
            "Available: ", shown, if (length(available) > 10) ", ...",
            call. = FALSE)
    }
    checkColumn(cells, outcomeCol, "outcomeCol")
    checkColumn(cells, "imageID", "imageID")

    # the relative-to-parent pipeline has no separate weighted/unweighted tables
    isRelative <- !"weighted" %in% names(results)
    dat <- results$proportions[results$proportions$cellType == cellType, ]
    yLab <- if (isRelative) {
        paste0("Proportion of ", unique(dat$group)[1], " cells")
    } else {
        "Proportion of all cells"
    }
    plotDat <- weightedLong(dat, "proportionUnweighted", "proportionWeighted",
        cells, outcomeCol)
    plotDat$type <- factor(plotDat$type, levels = c("Unweighted", "Weighted"))
    plotDat[[outcomeCol]] <- as.factor(plotDat[[outcomeCol]])
    plotDat <- plotDat[!is.na(plotDat$plotValue), , drop = FALSE]
    weightedBoxplot(plotDat, outcomeCol, yLab, cellType)
}

# The unweighted vs weighted boxplot shared by plotProportionBoxplot() and
# plotMetricBoxplot(): one box per type and outcome level, with one
# interactive point per image (hovering shows its imageID). `plotDat` is from
# weightedLong().
weightedBoxplot <- function(plotDat, outcomeCol, yLab, title) {
    ggplot(plotDat, aes(x = .data$type, y = .data$plotValue,
        colour = .data[[outcomeCol]])) +
        geom_boxplot(aes(group = interaction(.data[[outcomeCol]], .data$type)),
            outlier.shape = NA, width = 0.35, position = position_dodge(0.45)) +
        ggiraph::geom_point_interactive(aes(tooltip = .data$imageID,
            data_id = .data$imageID,
            group = interaction(.data[[outcomeCol]], .data$type)),
            position = position_jitterdodge(jitter.width = 0.1,
                dodge.width = 0.45),
            alpha = 0.6, size = 1.5) +
        labs(x = NULL, y = yLab, colour = NULL, title = title) +
        theme_classic()
}

# Co-localisation ----
# Plot of the results of runColocalisationModels().

#' Plot one pair's co-localisation by outcome
#'
#' Shows, for one pair of cell types or groups, the per-image co-localisation
#' measured with the original labels (unweighted) and with the bootstrap
#' confidence weights (weighted), for each level of the outcome. Each point is
#' one image.
#'
#' @param results Output of
#'   [runColocalisationModels()].
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   `imageID` and the outcome in `colData`.
#' @param fromType,toType The pair to show. An error is raised if the pair is
#'   not in the results; relative co-localisation (`"kontextual"`) only has
#'   pairs whose `toType` has a parent other than itself.
#' @param outcomeCol Name of the outcome column in `colData(spe)`.
#' @param metric The measure to plot: `"lfunction"` (co-localisation, the
#'   cross-type L-function) or `"kontextual"` (co-localisation relative to a
#'   parent population, the Kontextual measure).
#'
#' @return An interactive plot, a `girafe` widget that needs the ggiraph
#'   package.
#'
#' @examples
#' data(lungSubset)
#' data(lungColocalisation)
#' if (requireNamespace("ggiraph", quietly = TRUE) &&
#'     requireNamespace("gdtools", quietly = TRUE)) {
#'     plotMetricBoxplot(lungColocalisation, lungSubset, "Mast",
#'         "Proliferating FBs", outcomeCol = "disease_status")
#'     plotMetricBoxplot(lungColocalisation, lungSubset, "Mast",
#'         "Proliferating FBs", outcomeCol = "disease_status",
#'         metric = "kontextual")
#' }
#' @export
plotMetricBoxplot <- function(results, spe, fromType, toType, outcomeCol,
        metric = c("lfunction", "kontextual")) {
    checkInteractivePackages("plotMetricBoxplot")
    checkSpe(spe)
    metric <- match.arg(metric)
    p <- metricBoxplotGG(results, getCellData(spe), fromType, toType,
        outcomeCol, metric)
    interactivePlot(p)
}

# the ggplot behind plotMetricBoxplot(), before it is made interactive
metricBoxplotGG <- function(results, cells, fromType, toType, outcomeCol,
        metric) {
    relative <- metric == "kontextual"
    tables <- if (relative) {
        c("relativeUnweighted", "relativeWeighted")
    } else {
        c("unweighted", "weighted")
    }
    unweighted <- results[[tables[1]]]
    weighted <- results[[tables[2]]]
    needed <- c("imageID", "from", "to", "colocalisation")
    if (!is.data.frame(unweighted) || !is.data.frame(weighted) ||
        !all(needed %in% names(unweighted))) {
        stop("`results` must be the output of runColocalisationModels()",
            call. = FALSE)
    }
    checkPairType(fromType, "fromType")
    checkPairType(toType, "toType")
    checkColumn(cells, "imageID", "imageID")
    inPair <- function(d) d$from == fromType & d$to == toType
    if (!any(inPair(unweighted))) {
        hint <- if (relative) {
            paste0(" for metric = \"kontextual\". Relative co-localisation ",
                "only has pairs whose `toType` has a parent other than itself ",
                "(see ?runColocalisationModels)")
        }
        stop("The pair '", fromType, "' -> '", toType, "' is not in the ",
            "results", hint, call. = FALSE)
    }
    checkColumn(cells, outcomeCol, "outcomeCol")

    u <- unweighted[inPair(unweighted), c("imageID", "colocalisation")]
    w <- weighted[inPair(weighted), c("imageID", "colocalisation")]
    dat <- data.frame(imageID = u$imageID, value = u$colocalisation,
        weightedValue = w$colocalisation[match(u$imageID, w$imageID)])
    dat <- dat[dat$imageID %in% w$imageID, , drop = FALSE]
    plotDat <- weightedLong(dat, "value", "weightedValue", cells, outcomeCol)
    yLab <- if (relative) "Relative co-localisation" else "Co-localisation"
    weightedBoxplot(plotDat, outcomeCol, yLab,
        paste0(fromType, " \u2192 ", toType))
}

# Checks that `value`, the argument called `arg`, is a single cell type or
# group name.
checkPairType <- function(value, arg) {
    if (!is.character(value) || length(value) != 1 || is.na(value)) {
        stop("`", arg, "` must be a single cell type or group name",
            call. = FALSE)
    }
    invisible(TRUE)
}

# Niches ----
# Plots of the niches from runNicheIdentification() and testNicheOutcome().

# The relative enrichment of each cell type (rows) in each niche (columns): its
# share of the niche's cells divided by its share of all cells, so 1 means no
# enrichment.
nicheEnrichment <- function(cells, cellType, nicheCol) {
    tab <- table(cells[[cellType]], cells[[nicheCol]])
    tab / rowSums(tab) %*% t(colSums(tab)) * sum(tab)
}

#' Niche plots
#'
#' `plotNicheSpatial()` shows the cells of one niche in one image.
#' `nicheCompositionMap()` shows the relative enrichment of each cell type in
#' each niche as a bubble plot or heatmap. `nicheBarplot()` shows it for a
#' single niche. `nicheProportionBoxplot()` compares a niche's per-image
#' proportion between outcome groups.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   the cell type, image and coordinates of every cell in the form given by
#'   [getCellData()].
#' @param imageID The image to plot.
#' @param niche For `plotNicheSpatial()` and the proportion boxplot, the
#'   niche to show. For `nicheCompositionMap()` and `nicheBarplot()`, the
#'   niche assignments: either a vector of niches named by cell ID, such as
#'   `niches$nicheWeighted` from
#'   [runNicheIdentification()], or the name of
#'   a `colData` column holding them.
#' @param nicheCol The niche assignments for `plotNicheSpatial()`, in either
#'   of the forms taken by `niche` in `nicheCompositionMap()`.
#' @param cellSize Size of the highlighted cells in `plotNicheSpatial()`, the
#'   cells of the chosen niche. The grey background cells always have size 0.7.
#' @param cellType Column of `colData(spe)` with cell type labels.
#' @param type `"bubble"` or `"heatmap"`.
#' @param limit Lower and upper limit of the relative enrichment.
#' @param ... Further arguments passed to
#'   [pheatmap::pheatmap()].
#' @param target The niche to plot in `nicheBarplot()`.
#' @param proportions The `propWeighted` or `propUnweighted` element of
#'   the results of [testNicheOutcome()]: the table
#'   you pass decides whether the weighted or the unweighted niches are
#'   plotted.
#' @param outcomeCol Name of the per-image outcome column in `colData(spe)`.
#'
#' @return A `ggplot` object; for `nicheProportionBoxplot()` a `girafe`
#'   (interactive) widget; for `nicheCompositionMap(type = "heatmap")` the
#'   result of [pheatmap::pheatmap()].
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' set.seed(1)
#' niches <- runNicheIdentification(lungSubset, lungCellTypeTree, lungConfList,
#'     r = 50, k = 6)
#' plotNicheSpatial(lungSubset, "VUHD113", "niche1",
#'     nicheCol = niches$nicheWeighted)
#' nicheCompositionMap(lungSubset, niches$nicheWeighted)
#' nicheBarplot(lungSubset, niches$nicheWeighted, "niche1")
#' res <- testNicheOutcome(lungSubset, niches$nicheUnweighted,
#'     niches$nicheWeighted, outcomeCol = "disease_status", subject = "subject")
#' if (requireNamespace("ggiraph", quietly = TRUE) &&
#'     requireNamespace("gdtools", quietly = TRUE)) {
#'     nicheProportionBoxplot(res$propWeighted, "niche1", lungSubset,
#'         "disease_status")
#'     nicheProportionBoxplot(res$propUnweighted, "niche1", lungSubset,
#'         "disease_status")
#' }
#' @name nichePlots
#' @export
plotNicheSpatial <- function(spe, imageID, niche, nicheCol = "nicheWeighted",
        cellSize = 1) {
    checkCellSize(cellSize)
    checkSpe(spe)
    labelled <- withNicheColumn(getCellData(spe), nicheCol, "nicheCol")
    cells <- labelled$cells
    nicheCol <- labelled$col
    checkImage(cells, imageID)
    checkNiche(cells, niche, nicheCol,
        label = if (nicheCol == ".niche") "nicheCol" else nicheCol)
    cells <- cells[cells$imageID == imageID, , drop = FALSE]
    inNiche <- cells[[nicheCol]] %in% niche

    spatialPlot(cells[!inNiche, , drop = FALSE]) +
        geom_point(data = cells[inNiche, , drop = FALSE],
            aes(x = .data$x, y = .data$y, colour = .data$cellType),
            size = cellSize) +
        spatialFinish(paste0(imageID, " - ", niche), colour = "Cell type")
}

#' @rdname nichePlots
#' @export
nicheCompositionMap <- function(spe, niche, cellType = "cellType",
        type = c("bubble", "heatmap"), limit = c(0.33, 3), ...) {
    type <- match.arg(type)
    checkSpe(spe)
    labelled <- withNicheColumn(getCellData(spe), niche, "niche")
    checkColumn(labelled$cells, cellType, "cellType")
    tab <- nicheEnrichment(labelled$cells, cellType, labelled$col)
    tabClipped <- pmax(pmin(tab, limit[2]), limit[1])
    if (type == "heatmap") {
        return(pheatmap::pheatmap(tabClipped, cluster_cols = FALSE, ...))
    }

    # the bubble plot orders the cell types as the heatmap's clustering does
    ph <- pheatmap::pheatmap(tabClipped, cluster_cols = FALSE, silent = TRUE,
        ...)
    dat <- as.data.frame(tab)
    dat$cellType <- factor(dat$Var1, levels = levels(dat$Var1)[
        ph$tree_row$order])
    dat$niche <- dat$Var2
    dat$enrichment <- pmax(pmin(dat$Freq, limit[2]), limit[1])

    ggplot(dat, aes(x = .data$niche, y = .data$cellType,
        colour = .data$enrichment, size = .data$enrichment)) +
        geom_point() +
        scale_colour_gradient2(low = "#4575B4", mid = "grey90",
            high = "#D73027", midpoint = 1, guide = "legend") +
        theme_minimal() +
        labs(x = "Niche", y = "Cell type", colour = "Relative\nenrichment",
            size = "Relative\nenrichment") +
        theme(axis.text.x = element_text(angle = 90))
}

#' @rdname nichePlots
#' @export
nicheBarplot <- function(spe, niche, target, cellType = "cellType",
        limit = c(0.33, 3)) {
    checkSpe(spe)
    labelled <- withNicheColumn(getCellData(spe), niche, "niche")
    checkColumn(labelled$cells, cellType, "cellType")
    byColumn <- labelled$col != ".niche"
    tab <- nicheEnrichment(labelled$cells, cellType, labelled$col)
    if (!is.character(target) || length(target) != 1 ||
        !target %in% colnames(tab)) {
        niches <- paste(colnames(tab), collapse = ", ")
        stop("`target` must be one of the niches in `",
            if (byColumn) niche else "niche", "`: ", niches, call. = FALSE)
    }

    enrichment <- pmax(pmin(tab[, target], limit[2]), limit[1])
    dat <- data.frame(cellType = stats::reorder(rownames(tab), enrichment),
        enrichment = unname(enrichment))
    ggplot(dat, aes(x = .data$cellType, y = .data$enrichment)) +
        geom_col(fill = "#2a78d6", width = 0.7) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
        coord_flip() +
        labs(x = "Cell type", y = "Relative enrichment", title = if (byColumn) {
            paste0(target, " (", niche, ")")
        } else {
            target
        }) +
        theme_classic()
}

#' @rdname nichePlots
#' @export
nicheProportionBoxplot <- function(proportions, niche, spe, outcomeCol) {
    checkInteractivePackages("nicheProportionBoxplot")
    checkSpe(spe)
    p <- nicheProportionBoxplotGG(proportions, niche, getCellData(spe),
        outcomeCol)
    interactivePlot(p)
}

# the ggplot behind nicheProportionBoxplot(), before it is made interactive
nicheProportionBoxplotGG <- function(proportions, niche, cells, outcomeCol) {
    if (!is.data.frame(proportions) ||
        !all(c("imageID", "niche", "proportion") %in% names(proportions))) {
        stop("`proportions` must be the `propWeighted` or `propUnweighted` ",
            "element of the results of testNicheOutcome()", call. = FALSE)
    }
    if (!is.character(niche) || length(niche) != 1 ||
        !niche %in% proportions$niche) {
        niches <- paste(sort(unique(proportions$niche)), collapse = ", ")
        stop("`niche` must be one of the niches in `proportions`: ", niches,
            call. = FALSE)
    }
    checkColumn(cells, outcomeCol, "outcomeCol")
    checkColumn(cells, "imageID", "imageID")

    outcomeDat <- unique(cells[, c("imageID", outcomeCol)])
    plotDat <- proportions[proportions$niche == niche, , drop = FALSE]
    plotDat <- plotDat[plotDat$imageID %in% outcomeDat$imageID, , drop = FALSE]
    plotDat[[outcomeCol]] <- as.factor(outcomeDat[[outcomeCol]][
        match(plotDat$imageID, outcomeDat$imageID)])
    # testNicheOutcome() records whether the table is weighted or unweighted
    type <- attr(proportions, "nicheType")
    title <- paste0(niche, if (!is.null(type)) paste0(" (", type, ")"))

    ggplot(plotDat, aes(x = .data[[outcomeCol]], y = .data$proportion,
        colour = .data[[outcomeCol]])) +
        geom_boxplot(outlier.shape = NA, width = 0.5) +
        ggiraph::geom_point_interactive(aes(tooltip = .data$imageID,
            data_id = .data$imageID), position = position_jitter(width = 0.1),
            alpha = 0.7, size = 1.5) +
        labs(x = NULL, y = paste0(niche, " proportion"), colour = NULL,
            title = title) +
        theme_classic()
}
