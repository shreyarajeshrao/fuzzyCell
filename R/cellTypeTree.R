#' Cell type tree utilities
#'
#' Helpers for the cell type tree: a nested list describing how cell types are
#' organised into broad lineages and (optionally) groups of closely related
#' subtypes.
#' [bootstrapLineageConfidence()] refines
#' the cells within each group, and the analysis functions compare cell types,
#' groups and lineages.
#'
#' @details
#' A tree is a named list with one entry per broad lineage (these become the
#' `lineage` column). Each lineage is either a character vector of cell types,
#' or a list whose entries are:
#'
#' * a bare string: a cell type that is not in any group;
#' * a named `list(children = c(...))`: a **group** of closely
#'   related cell types, the unit that is bootstrap-refined. A group may contain
#'   only `children`;
#' * a named list of further entries: purely organisational nesting (a
#'   "wrapper") that is not a group itself, to any depth.
#'
#' A lineage can also itself be `list(children = c(...))`, in which case the
#' lineage name doubles as the group name. A cell type's group is its nearest
#' `list(children = ...)` ancestor, however many wrappers sit above it. In the
#' example below `CD8_Resting` is in the group `CD8+ T cells`, not
#' `T cells`.
#'
#' Every cell type must appear exactly once, and lineage, wrapper and group
#' names must be unique and different from every cell type name.
#'
#' ```
#' list(
#'     "Immune" = list(
#'         "T cells" = list(
#'             "CD8+ T cells" = list(children = c("CD8_Resting", "CD8_Memory")),
#'             "CD4+ T cells" = list(children = c("CD4_naive", "CD4_Memory"))),
#'         "NK",
#'         "Mast"),
#'     "Endothelial lineage" = c("Endothelial", "Activated endothelial"))
#' ```
#'
#' `resolveTreeLeaves()` returns every cell type below a node of the tree,
#' however deep. `autoParentOverrides()` returns the node directly above each
#' group. Both help choose the parent for relative co-localisation (see
#' `parentOverrides` in
#' [runColocalisationModels()]).
#'
#' @param cellTypeTree A nested list describing the cell type hierarchy (see
#'   Details and [lungCellTypeTree]).
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`), or
#'   a `data.frame`/`DataFrame` of per-cell data.
#' @param cellType Name of the column of `colData(spe)` (or of the data
#'   frame) holding each cell's cell type label.
#' @param nodeName A single string: the name of any node in `cellTypeTree`, a
#'   lineage, an organisational wrapper, a group or a cell type.
#'
#' @return
#' * `buildCellTypeMap()`: a list with `cellMap` (a `data.frame`
#'   with columns `cellType`, `group` and `lineage`, one row per cell
#'   type) and `groupsWithChildren` (names of the groups to be
#'   bootstrap-refined). An error is raised if the tree does not follow the
#'   rules in Details.
#' * `annotateCellTypeTree()`: `spe` with `group` and
#'   `lineage` columns added (a message reports any existing `group` or
#'   `lineage` column whose values are replaced). An error is raised if any
#'   cell's cell type is not in the tree (including `NA`).
#' * `resolveTreeLeaves()`: character vector of every cell type beneath
#'   `nodeName`.
#' * `autoParentOverrides()`: named character vector mapping each group
#'   that has a named ancestor to that ancestor, for use as `parentOverrides` in
#'   [runColocalisationModels()].
#'
#' @examples
#' data(lungCellTypeTree)
#' buildCellTypeMap(lungCellTypeTree)
#' resolveTreeLeaves(lungCellTypeTree, "Immune")
#' autoParentOverrides(lungCellTypeTree)
#'
#' data(lungSubset)
#' cells <- annotateCellTypeTree(getCellData(lungSubset), lungCellTypeTree)
#' table(cells$group, useNA = "ifany")
#' @name cellTypeTree
#' @export
buildCellTypeMap <- function(cellTypeTree) {
    if (!is.list(cellTypeTree) || length(cellTypeTree) == 0 ||
        is.null(names(cellTypeTree))) {
        stop("cellTypeTree must be a non-empty named list, with one named ",
            "entry per lineage (see ?cellTypeTree)", call. = FALSE)
    }
    unnamed <- which(is.na(names(cellTypeTree)) | !nzchar(names(cellTypeTree)))
    if (length(unnamed) > 0) {
        entries <- paste(unnamed, collapse = ", ")
        stop("Every top-level entry of cellTypeTree must be named (one per ",
            "lineage), but entry ", entries, " ",
            if (length(unnamed) == 1) "is" else "are",
            " unnamed. Give it a name, for example list(Other = \"NK\").",
            call. = FALSE)
    }
    invisible(Map(.validateTreeNode, cellTypeTree, names(cellTypeTree)))

    treeDF <- treeLeafTable(cellTypeTree)
    checkTreeNames(cellTypeTree, treeDF)
    list(
        cellMap = data.frame(cellType = treeDF$child, group = treeDF$group,
            lineage = treeDF$lineage),
        # the groups, the units bootstrapLineageConfidence() refines;
        # ungrouped cell types have group NA
        groupsWithChildren = unique(treeDF$group[!is.na(treeDF$group)])
    )
}

# One row per cell type of the tree: its lineage, its group (NA if it has
# none) and its name (`child`).
treeLeafTable <- function(cellTypeTree) {
    leafRows <- function(leaves, lineage, group) {
        data.frame(lineage = rep(lineage, length(leaves)),
            group = rep(group, length(leaves)), child = as.character(leaves))
    }
    # the rows of every cell type below `node`, a node of `lineage`
    walk <- function(node, lineage) {
        # a plain character vector: these cell types have no group
        if (is.character(node)) {
            return(leafRows(node, lineage, NA_character_))
        }
        # a lineage that is itself a group: its name is the group name
        if (.isGroup(node)) {
            return(leafRows(node[["children"]], lineage, lineage))
        }
        nms <- .entryNames(node)
        do.call(rbind, lapply(seq_along(node), function(i) {
            entry <- node[[i]]
            if (!nzchar(nms[i])) {
                # a bare cell type among its siblings, with no group
                leafRows(entry, lineage, NA_character_)
            } else if (.isGroup(entry)) {
                # a group: nms[i] is the group, its children the cell types
                leafRows(entry[["children"]], lineage, nms[i])
            } else {
                # purely organisational nesting (e.g. "T cells" wrapping
                # further groups): recurse without adding a group
                walk(entry, lineage)
            }
        }))
    }
    do.call(rbind, lapply(seq_along(cellTypeTree), function(i) {
        walk(cellTypeTree[[i]], lineage = names(cellTypeTree)[i])
    }))
}

# Checks the names in a tree: every cell type appears exactly once, lineage,
# wrapper and group names are unique, and none of them is also a cell type.
# `treeDF` is from treeLeafTable().
checkTreeNames <- function(cellTypeTree, treeDF) {
    # annotateCellTypeTree() matches cells to the tree by cell type, so a
    # repeated name would silently put all its cells in the first group
    repeated <- unique(treeDF$child[duplicated(treeDF$child)])
    if (length(repeated) > 0) {
        describe <- vapply(repeated, function(ct) {
            hits <- treeDF[treeDF$child == ct, ]
            where <- ifelse(is.na(hits$group),
                paste0("'", hits$lineage, "' (no group)"),
                paste0("'", hits$group, "'"))
            sprintf("'%s' (under %s)", ct, paste(where, collapse = " and "))
        }, character(1))
        details <- paste(describe, collapse = "; ")
        stop("Each cell type must appear exactly once in cellTypeTree, but ",
            details, if (length(repeated) == 1) " appears" else " appear",
            " more than once.", call. = FALSE)
    }

    # lineage, wrapper and group names must differ from cell type names: a
    # label that is both would always be read as the cell type, so the
    # lineage or group could never be analysed
    nodes <- .treeNodes(cellTypeTree)
    nodeNames <- c(nodes$name[is.na(nodes$parent)],
        nodes$name[!is.na(nodes$parent) & nodes$isList])
    clash <- intersect(nodeNames, treeDF$child)
    if (length(clash) > 0) {
        clashing <- paste0("'", clash, "'", collapse = ", ")
        stop(clashing, " ", if (length(clash) == 1) "is" else "are",
            " used both as a lineage/group name and as a cell type in ",
            "cellTypeTree. Names must be unique across the whole tree.",
            call. = FALSE)
    }
    repeatedNodes <- unique(nodeNames[duplicated(nodeNames)])
    if (length(repeatedNodes) > 0) {
        repeatedNames <- paste0("'", repeatedNodes, "'", collapse = ", ")
        stop("Lineage and group names must be unique in cellTypeTree, but ",
            repeatedNames, " ",
            if (length(repeatedNodes) == 1) "is" else "are",
            " used more than once.", call. = FALSE)
    }
    invisible(TRUE)
}

#' @rdname cellTypeTree
#' @export
annotateCellTypeTree <- function(spe, cellTypeTree, cellType = "cellType") {
    checkColumn(spe, cellType, "cellType")
    cellMap <- buildCellTypeMap(cellTypeTree)$cellMap

    # every cell must belong to a cell type in the tree: otherwise it would
    # quietly get NA for group and lineage and drop out of the analyses
    labels <- as.character(spe[[cellType]])
    notInTree <- setdiff(unique(labels), cellMap$cellType)
    if (length(notInTree) > 0) {
        nCells <- vapply(notInTree, function(ct) {
            if (is.na(ct)) {
                sum(is.na(labels))
            } else {
                sum(labels == ct, na.rm = TRUE)
            }
        }, integer(1))
        missingTypes <- paste0("'", notInTree, "' (", nCells, " cell",
            ifelse(nCells == 1, "", "s"), ")", collapse = ", ")
        stop(length(notInTree), " cell type", if (length(notInTree) > 1) "s",
            " (", sum(nCells), " cell", if (sum(nCells) != 1) "s",
            ") in the data ", if (length(notInTree) > 1) "are" else "is",
            " not in cellTypeTree: ", missingTypes,
            ". Add them to the tree, or remove these cells, before continuing.",
            call. = FALSE)
    }

    rows <- match(spe[[cellType]], cellMap$cellType)
    newValues <- list(group = cellMap$group[rows],
        lineage = cellMap$lineage[rows])

    # tell the user if a column of theirs is being replaced with different
    # values (re-annotating already annotated cells is silent)
    replaced <- names(newValues)[vapply(names(newValues), function(col) {
        col %in% cellDataNames(spe) &&
            !isTRUE(all.equal(as.character(spe[[col]]), newValues[[col]]))
    }, logical(1))]
    if (length(replaced) > 0) {
        columns <- paste0("'", replaced, "'", collapse = " and ")
        message("Replacing the existing ", columns, " column",
            if (length(replaced) > 1) "s", " in the data with ",
            if (length(replaced) > 1) "values" else "a column",
            " derived from cellTypeTree")
    }

    spe$group <- newValues$group
    spe$lineage <- newValues$lineage
    spe
}

# TRUE if `node` is a group: a list with an element named exactly "children"
# (exact matching, unlike `$`, which would also match e.g. "childrenX")
.isGroup <- function(node) {
    is.list(node) && "children" %in% names(node)
}

# The names of a node's entries, with "" for unnamed entries.
.entryNames <- function(node) {
    nms <- names(node)
    if (is.null(nms)) rep("", length(node)) else nms
}

# Checks that the part of the tree below one node is well formed, and stops
# with a message saying where if not. `path` is the vector of names leading to
# `node`, used only for the message.
.validateTreeNode <- function(node, path) {
    where <- paste(path, collapse = " > ")
    if (is.character(node)) {
        if (length(node) == 0 || anyNA(node) || any(!nzchar(node))) {
            stop("'", where, "' has a missing or empty cell type name; every ",
                "cell type must be a non-empty string", call. = FALSE)
        }
        return(invisible())
    }
    if (!is.list(node)) {
        what <- if (is.null(node)) "NULL" else paste("a", class(node)[1])
        stop("'", where, "' must be a character vector of cell types or a ",
            "list, not ", what, call. = FALSE)
    }

    nms <- .entryNames(node)
    if (.isGroup(node)) {
        .validateGroup(node, nms, where)
        return(invisible())
    }
    # a loop, not lapply(), so that the first problem stops the check
    for (i in seq_along(node)) {
        if (is.list(node[[i]]) && !nzchar(nms[i])) {
            stop("Entry ", i, " of '", where, "' is an unnamed list; nested ",
                "entries need a name", call. = FALSE)
        }
        entryName <- if (nzchar(nms[i])) nms[i] else sprintf("[[%d]]", i)
        .validateTreeNode(node[[i]], c(path, entryName))
    }
    invisible()
}

# Checks that a group holds only one `children` entry: a non-empty character
# vector of cell type names. `where` names the group in the message.
.validateGroup <- function(node, nms, where) {
    extra <- setdiff(seq_along(node), which(nms == "children"))
    if (length(extra) > 0 || sum(nms == "children") > 1) {
        others <- if (length(extra) > 0) {
            labels <- ifelse(nzchar(nms[extra]), nms[extra],
                sprintf("[[%d]]", extra))
            sprintf(", but it also has %d other entr%s (%s)", length(extra),
                if (length(extra) == 1) "y" else "ies",
                paste0("'", labels, "'", collapse = ", "))
        } else {
            ""
        }
        stop("The group '", where, "' may only contain `children = c(...)`",
            others, ". Put other cell types next to the group, not inside it.",
            call. = FALSE)
    }
    children <- node[["children"]]
    if (!is.character(children) || length(children) == 0 ||
        anyNA(children) || any(!nzchar(children))) {
        stop("`children` of '", where, "' must be a non-empty character ",
            "vector of cell type names, without missing or empty names",
            call. = FALSE)
    }
    invisible()
}

# Every cell type below a node of the tree, however deep.
.treeLeaves <- function(node) {
    if (is.character(node)) return(node)
    if (.isGroup(node)) return(node[["children"]])
    unlist(lapply(node, .treeLeaves), use.names = FALSE)
}

# One row per named node of a (valid) tree, in depth-first order: lineages,
# wrappers, groups and named vectors of cell types. Columns: `name`, `parent`
# (the node it sits in; NA for a lineage), `isList` (FALSE for a named vector
# of cell types) and `leaves` (a list column of the cell types below it).
.treeNodes <- function(cellTypeTree) {
    visit <- function(node, name, parent) {
        here <- list(list(name = name, parent = parent, isList = is.list(node),
            leaves = .treeLeaves(node)))
        if (!is.list(node) || .isGroup(node)) return(here)
        nms <- .entryNames(node)
        below <- lapply(which(nzchar(nms)), function(i) {
            visit(node[[i]], nms[i], name)
        })
        c(here, unlist(below, recursive = FALSE))
    }
    rows <- unlist(lapply(seq_along(cellTypeTree), function(i) {
        visit(cellTypeTree[[i]], names(cellTypeTree)[i], NA_character_)
    }), recursive = FALSE)
    nodes <- data.frame(
        name = vapply(rows, `[[`, "", "name"),
        parent = vapply(rows, `[[`, "", "parent"),
        isList = vapply(rows, `[[`, TRUE, "isList"))
    nodes$leaves <- lapply(rows, `[[`, "leaves")
    nodes
}

#' @rdname cellTypeTree
#' @export
resolveTreeLeaves <- function(cellTypeTree, nodeName) {
    if (!is.character(nodeName) || length(nodeName) != 1 || is.na(nodeName)) {
        stop("`nodeName` must be a single, non-missing string", call. = FALSE)
    }
    allLeaves <- buildCellTypeMap(cellTypeTree)$cellMap$cellType
    if (nodeName %in% allLeaves) {
        return(nodeName)
    }
    nodes <- .treeNodes(cellTypeTree)
    hit <- match(nodeName, nodes$name)
    if (is.na(hit)) {
        stop("No node named '", nodeName, "' found in cellTypeTree",
            call. = FALSE)
    }
    unique(nodes$leaves[[hit]])
}

#' @rdname cellTypeTree
#' @export
autoParentOverrides <- function(cellTypeTree) {
    groups <- buildCellTypeMap(cellTypeTree)$groupsWithChildren
    nodes <- .treeNodes(cellTypeTree)
    # a lineage that is itself a group has no parent (NA)
    parents <- stats::setNames(nodes$parent[match(groups, nodes$name)], groups)
    parents[!is.na(parents)]
}
