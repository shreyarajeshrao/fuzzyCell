# Helpers shared by several files of the package. Helpers used by one file
# only are kept in that file.

#' Per-cell data frame from a SummarizedExperiment
#'
#' Returns `colData(spe)` as a `data.frame`. For a
#' `SpatialExperiment` whose coordinates are stored in `spatialCoords()`
#' rather than in `colData()`, the coordinates are added as columns `x`
#' and `y`. The spatial functions in fuzzyCell expect `cellID`,
#' `imageID`, `cellType`, `x` and `y` columns.
#'
#' @param spe A `SummarizedExperiment` or `SpatialExperiment`.
#'
#' @return A `data.frame` with one row per cell.
#'
#' @examples
#' data(lungSubset)
#' head(getCellData(lungSubset))
#' @export
getCellData <- function(spe) {
    cells <- as.data.frame(SummarizedExperiment::colData(spe))
    if (methods::is(spe, "SpatialExperiment") &&
        !all(c("x", "y") %in% names(cells))) {
        xy <- SpatialExperiment::spatialCoords(spe)
        cells$x <- xy[, 1]
        cells$y <- xy[, 2]
    }
    cells
}

# A copy of the per-cell data with the standard columns the analysis code reads
# (cellID, cellType, imageID, x, y) filled from the columns the user named.
# `cols` is a named list: names are the standard columns, values the user's
# column names (each checked to be a single string). Checks that every column
# exists, that the coordinates are numeric and complete, and that cell IDs are
# unique.
standardiseCells <- function(cells, cols) {
    # a loop, so that the first missing column stops the check
    for (std in names(cols)) checkColumn(cells, cols[[std]], std)
    cells[names(cols)] <- cells[unlist(cols, use.names = FALSE)]

    for (coord in intersect(c("x", "y"), names(cols))) {
        if (!is.numeric(cells[[coord]]) || anyNA(cells[[coord]])) {
            stop("`", coord, " = \"", cols[[coord]], "\"` must be a numeric ",
                "column without missing values", call. = FALSE)
        }
    }
    if ("cellID" %in% names(cols) &&
        (anyNA(cells$cellID) || anyDuplicated(cells$cellID))) {
        example <- if (anyDuplicated(cells$cellID)) {
            sprintf(" (e.g. '%s' is repeated)",
                cells$cellID[anyDuplicated(cells$cellID)])
        }
        stop("The `cellID` column '", cols[["cellID"]], "' must contain ",
            "unique, non-missing cell IDs", example, call. = FALSE)
    }
    cells
}

# Per-cell data without any group / lineage columns from an earlier annotation,
# so that annotating it again does not report a replaced column.
withoutAnnotation <- function(cells) {
    cells[, !names(cells) %in% c("group", "lineage"), drop = FALSE]
}

# Names of the per-cell columns of a SummarizedExperiment (its colData) or of a
# data.frame / DataFrame.
cellDataNames <- function(x) {
    if (methods::is(x, "SummarizedExperiment")) {
        names(SummarizedExperiment::colData(x))
    } else {
        names(x)
    }
}

# Argument checks ----

# Checks that `spe` is a SummarizedExperiment (which includes SpatialExperiment
# and SingleCellExperiment).
checkSpe <- function(spe) {
    if (!methods::is(spe, "SummarizedExperiment")) {
        stop("`spe` must be a SummarizedExperiment, for example a ",
            "SpatialExperiment", call. = FALSE)
    }
    invisible(spe)
}

# Checks that `value`, the argument called `arg`, is a single string naming a
# column of x (colData of a SummarizedExperiment, or a data.frame).
checkColumn <- function(x, value, arg) {
    where <- if (methods::is(x, "SummarizedExperiment")) {
        "colData(spe)"
    } else {
        "the data"
    }
    if (!is.character(value) || length(value) != 1 || is.na(value)) {
        stop("`", arg, "` must be a single string naming a column of ", where,
            call. = FALSE)
    }
    available <- cellDataNames(x)
    if (!value %in% available) {
        shown <- paste(utils::head(available, 10), collapse = ", ")
        stop("`", arg, " = \"", value, "\"` is not a column of ", where,
            ". Available columns: ", shown,
            if (length(available) > 10) ", ...", call. = FALSE)
    }
    invisible(value)
}

# Checks that `value`, the argument called `arg`, is a single string naming an
# assay of the SummarizedExperiment `spe`.
checkAssay <- function(spe, value, arg = "assayName") {
    if (!is.character(value) || length(value) != 1 || is.na(value)) {
        stop("`", arg, "` must be a single string naming an assay of spe",
            call. = FALSE)
    }
    available <- SummarizedExperiment::assayNames(spe)
    if (!value %in% available) {
        shown <- paste(available, collapse = ", ")
        stop("`", arg, " = \"", value, "\"` is not an assay of spe. ",
            "Available assays: ", shown, call. = FALSE)
    }
    invisible(value)
}

# Checks that `cores` is a single whole number of at least 1.
checkCores <- function(cores, arg = "cores") {
    if (!is.numeric(cores) || length(cores) != 1 || is.na(cores) ||
        cores < 1 || cores != round(cores)) {
        stop("`", arg, "` must be a single whole number of at least 1",
            call. = FALSE)
    }
    invisible(cores)
}

# Checks that `r` is a non-empty vector of positive, finite radii.
checkRadii <- function(r, arg = "r") {
    if (!is.numeric(r) || length(r) == 0 || anyNA(r) || any(!is.finite(r)) ||
        any(r <= 0)) {
        stop("`", arg, "` must be a vector of one or more positive numbers ",
            "(radii)", call. = FALSE)
    }
    invisible(r)
}

# Outcome and subjects ----

# Checks the per-image outcome that the analyses test against: the column
# exists, has a single value per image, and has exactly two levels. If
# `subject` is given, also checks that every image belongs to one subject;
# without it, says that images are treated as independent.
checkOutcome <- function(cells, outcomeCol, subject = NULL) {
    checkColumn(cells, outcomeCol, "outcomeCol")
    checkColumn(cells, "imageID", "imageID")
    if (!is.null(subject)) checkSubjects(cells, subject)
    outcome <- cells[[outcomeCol]]

    # one value per image: each image is one observation
    perImage <- split(outcome, cells$imageID)
    valuesPerImage <- lapply(perImage, function(v) unique(v[!is.na(v)]))
    mixed <- names(valuesPerImage)[lengths(valuesPerImage) > 1]
    if (length(mixed) > 0) {
        examples <- paste(vapply(utils::head(mixed, 5), function(id) {
            sprintf("%s (%s)", id, paste(valuesPerImage[[id]], collapse = ", "))
        }, character(1)), collapse = "; ")
        stop("`outcomeCol = \"", outcomeCol, "\"` must have a single value ",
            "per image, but ", length(mixed), " image",
            if (length(mixed) > 1) "s", " ",
            if (length(mixed) > 1) "have" else "has", " more than one: ",
            examples, if (length(mixed) > 5) "; ...",
            ". The models treat each image as one observation, so please ",
            "make sure there is exactly one outcome value per image. Several ",
            "images from one subject are handled with the `subject` argument.",
            call. = FALSE)
    }

    # exactly two conditions
    levels <- sort(unique(outcome[!is.na(outcome)]))
    shownLevels <- paste(utils::head(levels, 6), collapse = ", ")
    if (length(levels) > 2) {
        stop("`outcomeCol = \"", outcomeCol, "\"` has ", length(levels),
            " levels (", shownLevels, if (length(levels) > 6) ", ...", "). ",
            "The models currently compare two conditions only. Please subset ",
            "the data to two conditions before running.", call. = FALSE)
    }
    if (length(levels) < 2) {
        only <- if (length(levels) == 1) sprintf(" (only %s)", shownLevels)
        stop("`outcomeCol = \"", outcomeCol, "\"` has fewer than two levels",
            only, ". Two conditions are needed to compare.", call. = FALSE)
    }
    messageIfNoSubject(subject)
    invisible(levels)
}

# Checks the `subject` column: no missing values, and every image belongs to a
# single subject. A subject may have several images; the analyses account for
# them (spicyR's and lmerTest's mixed models, or pooling).
checkSubjects <- function(cells, subject) {
    checkColumn(cells, subject, "subject")
    checkColumn(cells, "imageID", "imageID")
    ids <- as.character(cells[[subject]])
    if (anyNA(ids)) {
        noSubject <- unique(as.character(cells$imageID)[is.na(ids)])
        shown <- paste(utils::head(noSubject, 5), collapse = ", ")
        stop("The `subject` column '", subject, "' has missing values (in ",
            length(noSubject), " image", if (length(noSubject) > 1) "s", ": ",
            shown, if (length(noSubject) > 5) ", ...",
            "). Every cell needs a subject so that images from the same ",
            "subject can be detected.", call. = FALSE)
    }

    map <- unique(data.frame(imageID = as.character(cells$imageID),
        subject = ids))
    # every image must belong to one subject
    subjectsPerImage <- split(map$subject, map$imageID)
    shared <- subjectsPerImage[lengths(subjectsPerImage) > 1]
    if (length(shared) > 0) {
        examples <- paste(vapply(utils::head(names(shared), 5), function(img) {
            sprintf("%s (%s)", img,
                paste(utils::head(shared[[img]], 4), collapse = ", "))
        }, character(1)), collapse = "; ")
        stop("Every image must belong to a single subject, but ",
            length(shared), " image", if (length(shared) > 1) "s", " ",
            if (length(shared) > 1) "have" else "has",
            " more than one value in the `subject` column '", subject, "': ",
            examples, if (length(shared) > 5) "; ...", ".", call. = FALSE)
    }
    invisible(TRUE)
}

# One-line note that, without a `subject` column, images from one subject
# cannot be recognised.
messageIfNoSubject <- function(subject) {
    if (is.null(subject)) {
        message("`subject` was not given, so every image is treated as an ",
            "independent sample (from a different subject).")
    }
    invisible()
}

# TRUE if `subject` is given and some subject has several images.
hasRepeatedImages <- function(cells, subject) {
    if (is.null(subject)) return(FALSE)
    perImage <- unique(data.frame(imageID = as.character(cells$imageID),
        subject = as.character(cells[[subject]])))
    anyDuplicated(perImage$subject) > 0
}

# Stops with an install hint if lmerTest, needed for the mixed models when
# subjects have several images, is not installed.
needLmerTest <- function() {
    if (!requireNamespace("lmerTest", quietly = TRUE)) {
        stop("Some subjects have several images, which needs a mixed ",
            "model from the lmerTest package. Install it with ",
            "BiocManager::install(\"lmerTest\").", call. = FALSE)
    }
    invisible(TRUE)
}

# A linear mixed model of the column `response` of `dat` on the outcome, with
# a random intercept per subject (`.subject`), fitted with lmerTest::lmer() as
# spicyR does. Returns the estimate for the second level of the outcome, its t
# statistic, its Satterthwaite p-value and whether the fit was singular (1) or
# not (0), or NULL if the model cannot be fitted. lmerTest's "boundary
# (singular) fit" message is caught and recorded, not printed: the callers
# summarise it once (see reportSingular()). Other messages and warnings, such
# as convergence problems, reach the user.
mixedModelTest <- function(dat, response, outcomeCol) {
    formula <- stats::as.formula(paste0("`", response, "` ~ `", outcomeCol,
        "` + (1 | .subject)"))
    state <- new.env()
    state$singular <- FALSE
    fit <- tryCatch(withCallingHandlers(lmerTest::lmer(formula, data = dat),
        message = function(m) {
            if (grepl("singular", conditionMessage(m), fixed = TRUE)) {
                state$singular <- TRUE
                invokeRestart("muffleMessage")
            }
        }), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    coefs <- summary(fit)$coefficients
    if (nrow(coefs) < 2) return(NULL)
    c(estimate = coefs[2, "Estimate"], statistic = coefs[2, "t value"],
        p.value = coefs[2, "Pr(>|t|)"], singular = as.numeric(state$singular))
}

# Prints one message saying how many of the mixed models in the table `out`
# (column `singular`, NA where no mixed model was fitted) had a singular fit,
# and returns the table without that column. `what` names the tests, for
# example "unweighted cell type proportions".
reportSingular <- function(out, what) {
    fitted <- sum(!is.na(out$singular))
    nSingular <- sum(out$singular == 1, na.rm = TRUE)
    if (nSingular > 0) {
        message("Mixed models for the ", what, ": for ", nSingular,
            " of the ", fitted, " models the between-subject variance was ",
            "estimated as zero ",
            "(a singular fit), so those models are equivalent to linear ",
            "models.")
    }
    out$singular <- NULL
    out
}

# Classifier helpers, shared by bootstrapping and cross-fitting ----

# Log1p-normalises and z-scores every row (cell) of a cells x genes count
# matrix in one vectorised pass. Cells with zero total count come out
# all-zero.
scaleCountMat <- function(mat) {
    totals <- rowSums(mat)
    norm <- log1p(mat / pmax(totals, 1) * 1e4)
    norm[totals == 0, ] <- 0

    means <- rowMeans(norm)
    centred <- norm - means
    sds <- sqrt(rowSums(centred^2) / (ncol(mat) - 1))
    sds[sds == 0] <- 1 # avoid dividing by zero for zero-variance rows
    centred / sds
}

# TRUE for the columns of a matrix whose values are not all the same.
varyingColumns <- function(mat) {
    if (nrow(mat) < 2) return(rep(FALSE, ncol(mat)))
    colSums(mat != rep(mat[1, ], each = nrow(mat))) > 0
}

# Checks a cells x genes count matrix and its labels before classification: one
# label per row, and non-negative counts with no missing values. The resampling
# treats each cell's total as a number of transcripts, so unless
# `correctedCounts` is TRUE the counts must also be whole numbers (within a
# small tolerance for whole numbers stored as doubles).
checkCountMatrix <- function(countMat, labels, correctedCounts = FALSE) {
    if (nrow(countMat) != length(labels)) {
        stop("countMat has ", nrow(countMat), " rows but labels has ",
            length(labels), " elements; they must match (one row per cell)",
            call. = FALSE)
    }
    invalid <- anyNA(countMat) || any(countMat < 0)
    if (correctedCounts && invalid) {
        stop("countMat must contain non-negative counts, with no missing ",
            "values", call. = FALSE)
    }
    if (!correctedCounts &&
        (invalid || any(abs(countMat - round(countMat)) > 1e-6))) {
        stop("countMat must contain raw non-negative integer counts, with no ",
            "missing values. Did you pass normalised or log-transformed data? ",
            "Counts corrected for ambient RNA, which need not be whole ",
            "numbers, are accepted with correctedCounts = TRUE.", call. = FALSE)
    }
    invisible(TRUE)
}

# Checks that every class has at least 2 cells; `reason` says why, and follows
# "Each class needs at least 2 cells" in the error.
checkClassSizes <- function(labels, reason) {
    sizes <- table(labels)
    tooSmall <- names(sizes)[sizes < 2]
    if (length(tooSmall) > 0) {
        details <- paste0("'", tooSmall, "' has ", sizes[tooSmall],
            collapse = ", ")
        stop("Each class needs at least 2 cells", reason, ", but ", details,
            ". Drop or merge these cell types.", call. = FALSE)
    }
    invisible(sizes)
}

# The count matrix with its genes named by position if they have no names: the
# logistic fits find their genes by name.
withGeneNames <- function(countMat) {
    if (is.null(colnames(countMat))) {
        colnames(countMat) <- paste0("gene", seq_len(ncol(countMat)))
    }
    countMat
}

# The scaled counts (scaleCountMat()) and the genes that vary across the cells,
# which are the ones a classifier can use. Stops if no gene varies.
classifierFeatures <- function(countMat) {
    scaled <- scaleCountMat(countMat)
    keep <- varyingColumns(scaled)
    keep <- keep & !is.na(keep)
    if (!any(keep)) {
        stop("Every gene has zero variance across these cells, so there is ",
            "nothing to classify on. Check that countMat is a cells x genes ",
            "matrix of counts.", call. = FALSE)
    }
    list(scaled = scaled, keep = keep)
}

# A training sample of the cells `idx` with the same number of cells of each
# class: as many as the smallest class has, or `nPerClass` if that is more, in
# which case the smaller classes are drawn with replacement.
balancedSample <- function(idx, labels, nPerClass) {
    sizes <- table(labels[idx])
    targetSize <- max(min(sizes), nPerClass)
    unlist(lapply(unique(labels[idx]), function(lab) {
        pool <- idx[labels[idx] == lab]
        # sample.int() draws positions in the pool: sample() on a single
        # value k would draw from 1:k instead
        pool[sample.int(length(pool), size = targetSize,
            replace = length(pool) < targetSize)]
    }))
}
