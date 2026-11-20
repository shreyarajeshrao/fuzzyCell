#' Cell type weights from bootstrap confidences
#'
#' Builds the soft cell type membership of every cell, the weights that the
#' weighted analyses pass to spicyR and lisaClust: one column per leaf cell type
#' and per group of the tree. A cell of a bootstrapped group has, for each
#' subtype of the group, its bootstrap confidence in that subtype; every other
#' cell has weight 1 for its own cell type. A group's column is 1 for each of
#' its cells, since the bootstrap reclassifies cells only within their group.
#'
#' @inheritParams runColocalisationModels
#' @param weighted `TRUE` for the bootstrap confidences, `FALSE` for the
#'   original labels (0/1 weights) on the same columns.
#' @param cellType,cellID Names of the `colData(spe)` columns holding each
#'   cell's cell type label and a unique cell ID (matching the `cellID` column
#'   of `confList`).
#'
#' @return A numeric matrix with one row per cell of `spe` (in its order) and
#'   one named column per leaf cell type and group.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' w <- confidenceWeights(lungSubset, lungCellTypeTree, lungConfList)
#' head(w[, c("Arteriole", "Venous", "Endothelial lineage")])
#' @export
confidenceWeights <- function(spe, cellTypeTree, confList, weighted = TRUE,
        cellType = "cellType", cellID = "cellID") {
    checkSpe(spe)
    cells <- standardiseCells(getCellData(spe), list(cellID = cellID,
        cellType = cellType))
    weightMatrix(cells, cellTypeTree, confList, weighted)
}

# The cells x labels weight matrix behind confidenceWeights(), also used by
# runColocalisationModels() and runNicheIdentification(). `cells` is a per-cell
# table with cellID and cellType columns. Each cell has weight 1 for its own
# cell type and its group; with `weighted = TRUE`, a cell of a bootstrapped
# group has its confidences in place of 0/1 for that group's subtypes.
weightMatrix <- function(cells, cellTypeTree, confList, weighted) {
    # stops if a cell's type is not in the tree
    groupOf <- annotateCellTypeTree(withoutAnnotation(cells),
        cellTypeTree)$group
    treeInfo <- buildCellTypeMap(cellTypeTree)
    leaves <- treeInfo$cellMap$cellType
    groups <- treeInfo$groupsWithChildren
    labels <- c(leaves, setdiff(groups, leaves))
    w <- matrix(0, nrow(cells), length(labels), dimnames = list(NULL, labels))
    cellRows <- seq_len(nrow(cells))
    w[cbind(cellRows, match(as.character(cells$cellType), labels))] <- 1
    inGroup <- !is.na(groupOf)
    w[cbind(cellRows[inGroup], match(groupOf[inGroup], labels))] <- 1
    if (!weighted) return(w)

    # one group at a time: each fills its own cells' rows and its own
    # subtypes' columns, from its own confidence table
    for (g in intersect(names(confList), groups)) {
        conf <- confList[[g]]
        idx <- which(groupOf == g)
        rows <- match(as.character(cells$cellID[idx]),
            as.character(conf$cellID))
        if (anyNA(rows)) {
            stop("confList[[\"", g, "\"]] has no confidence for ",
                sum(is.na(rows)), " of the group's cells (matched by cell ID)",
                call. = FALSE)
        }
        # clear the group's 0/1 subtype labels, so only the confidences remain
        w[idx, treeInfo$cellMap$cellType[treeInfo$cellMap$group %in% g]] <- 0
        subtypes <- intersect(setdiff(names(conf), "cellID"), leaves)
        w[idx, subtypes] <- as.matrix(conf[rows, subtypes, drop = FALSE])
    }
    w
}
