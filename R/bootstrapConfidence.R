# Trains a classifier on the scaled counts `x` (cells x genes) of the training
# cells, with labels `y`. `ntree` is only used by the random forest, which is
# grown as 25 sub-forests on the workers of BPPARAM and merged with
# randomForest::combine(). The number of sub-forests does not depend on the
# number of workers, so neither does the forest.
trainClassifier <- function(x, y, classifier, ntree = 500,
        BPPARAM = parallelParam(1), nPC = 50) {
    y <- factor(y)
    if (classifier %in% c("logistic", "pcaLogistic")) {
        # pcaLogistic fits on the training cells' first nPC principal
        # components in place of the genes
        proj <- NULL
        k <- min(nPC, ncol(x) - 1, nrow(x) - 1)
        # glmnet needs at least two predictors, so with fewer than three genes
        # the genes are used as they are
        if (classifier == "pcaLogistic" && k >= 2) {
            pc <- principalComponents(x, k)
            proj <- list(genes = colnames(x), center = pc$center,
                rotation = pc$rotation)
            x <- pc$x
            colnames(x) <- paste0("PC", seq_len(ncol(x)))
        }

        cv <- glmnet::cv.glmnet(x, y, family = "multinomial", alpha = 0,
            nfolds = 3, type.measure = "deviance")
        co <- stats::coef(cv, s = "lambda.min")
        beta <- do.call(cbind, lapply(co, function(b) as.numeric(b[-1])))
        colnames(beta) <- names(co)
        return(list(a0 = vapply(co, function(b) b[1], 0), beta = beta,
            genes = colnames(x), proj = proj))
    }
    if (classifier == "rf") {
        # split ntree into 25 sub-forests of nearly equal size
        sizes <- diff(round(seq(0, ntree, length.out = min(25, ntree) + 1)))
        forests <- parLapply(sizes, function(nt) {
            randomForest::randomForest(x, y, ntree = nt)
        }, BPPARAM)
        return(do.call(randomForest::combine, forests))
    }
    switch(classifier,
        lda = MASS::lda(x, y),
        svm = e1071::svm(x, y, kernel = "radial")
    )
}

# The first k principal components of x (centred, not scaled), in the format of
# irlba::prcomp_irlba(). For up to 5,000 genes they are computed exactly, from
# an eigendecomposition of the genes' cross-product matrix. This is also faster
# than irlba: on 10,000 cells it took 0.8 s instead of 11 s with 1,000 genes,
# and 34 s instead of 84 s with 5,000 genes. Beyond 5,000 genes the
# eigendecomposition becomes too slow (its cost grows with the cube of the
# number of genes), so irlba is used instead.
principalComponents <- function(x, k) {
    if (ncol(x) > 5000) {
        return(irlba::prcomp_irlba(x, n = k, center = TRUE, scale. = FALSE))
    }
    center <- colMeans(x)
    x <- sweep(x, 2, center)
    rotation <- eigen(crossprod(x), symmetric = TRUE)$vectors[, seq_len(k),
        drop = FALSE]
    list(center = center, rotation = rotation, x = x %*% rotation)
}

# Precomputes what resampleHits() needs to resample the counts of every cell.
# Takes a cells x genes count matrix and the block size, in transcripts.
# Returns a list with:
# - nz: the row and column of every nonzero count, ordered by cell, then gene;
# - ends: the running total of those counts;
# - first, last: each cell's first and last row in nz;
# - start, totals: where each cell's counts start in `ends`, and their total;
# - draws: the number of transcripts each cell's resample draws (its total,
#   rounded);
# - blocks: the cells grouped into blocks of about `blockSize` transcripts,
#   resampled one block at a time to bound memory;
# - dim, dimnames: those of countMat.
resampleLayout <- function(countMat, blockSize = 2^22) {
    nz <- which(countMat > 0, arr.ind = TRUE)
    nz <- nz[order(nz[, 1], nz[, 2]), , drop = FALSE]
    ends <- cumsum(countMat[nz])
    # for each cell: its first and last row in nz, and where its counts start
    # in the running total `ends`
    last <- cumsum(tabulate(nz[, 1], nrow(countMat)))
    first <- c(1L, last[-length(last)] + 1L)
    start <- c(0, ends)[first]
    totals <- c(0, ends)[last + 1] - start
    draws <- round(totals)
    list(nz = nz, ends = ends, first = first, last = last, start = start,
        totals = totals, draws = draws,
        blocks = split(seq_along(draws), cumsum(draws) %/% blockSize),
        dim = dim(countMat), dimnames = dimnames(countMat))
}

# One resample: the new count of each nonzero entry in layout$nz (genes with
# a zero count stay zero). A cell's resample draws as many transcripts as its
# total count (rounded), each from gene j with probability count_j / total;
# counts that are not whole numbers are used as they are. The draws do not
# depend on how the cells are split into blocks.
resampleHits <- function(layout) {
    hits <- integer(nrow(layout$nz))
    for (cells in layout$blocks) {
        cell <- rep.int(cells, layout$draws[cells])
        if (length(cell) == 0) next
        # each transcript is a uniform random position along the cell's
        # counts, assigned to the gene whose stretch it falls in; clamping to
        # the cell's own entries guards against rounding at the ends
        pos <- layout$start[cell] +
            stats::runif(length(cell)) * layout$totals[cell]
        entry <- findInterval(pos, layout$ends) + 1
        entry <- pmin(pmax(entry, layout$first[cell]), layout$last[cell])
        # a block's cells cover consecutive rows of nz, so only that range of
        # `hits` needs filling
        lo <- min(entry)
        hits[lo:max(entry)] <- tabulate(entry - lo + 1L, max(entry) - lo + 1L)
    }
    hits
}

# The resampled counts as a cells x genes matrix.
resampleCounts <- function(layout) {
    out <- matrix(0L, layout$dim[1], layout$dim[2],
        dimnames = layout$dimnames)
    out[layout$nz] <- resampleHits(layout)
    out
}

# Turns a logistic fit into a single linear map from the scaled counts of all
# genes (scaleCountMat()) to the linear predictors: eta = z %*% W + offset.
# Genes the fit does not use get zero weight. For pcaLogistic, the projection
# onto the principal components is folded into W, and the centring into the
# offset.
linearMap <- function(fit, keep) {
    # beta has one row per fit$genes, and the rotation one row per
    # fit$proj$genes, in the same order
    W <- fit$beta
    offset <- fit$a0
    genes <- fit$genes
    if (!is.null(fit$proj)) {
        W <- fit$proj$rotation %*% W
        offset <- offset - drop(fit$proj$center %*% W)
        genes <- fit$proj$genes
    }
    full <- matrix(0, length(keep), ncol(W),
        dimnames = list(NULL, colnames(fit$beta)))
    full[match(genes, names(keep)), ] <- W
    list(W = full, offset = offset, colSumW = colSums(full))
}

# The linear predictors of a resample, computed from its nonzero entries
# only. They equal (up to rounding) those from scaling the full resample
# with scaleCountMat(). A cell's scaled value is z = (l - m) / s, where
# l = log1p(count / total * 1e4) is zero wherever the count is, and m and s
# are the mean and SD of l over all genes. So
# z %*% W = (l %*% W - m * colSums(W)) / s,
# which only needs the nonzero values of l.
resampleEta <- function(layout, hits, map) {
    nGenes <- layout$dim[2]
    cell <- layout$nz[, 1]
    l <- log1p(hits / pmax(layout$draws, 1)[cell] * 1e4)
    # genes x cells, column-compressed: entries ordered by cell, then gene
    lt <- Matrix::sparseMatrix(i = layout$nz[, 2], p = c(0L, layout$last),
        x = l, dims = c(nGenes, layout$dim[1]))
    eta <- as.matrix(Matrix::crossprod(lt, map$W))
    m <- Matrix::colSums(lt) / nGenes
    s <- sqrt(pmax(Matrix::colSums(lt^2) - nGenes * m^2, 0) /
        (nGenes - 1))
    s[s == 0] <- 1
    eta <- (eta - outer(m, map$colSumW)) / s
    sweep(eta, 2, map$offset, "+")
}

# The linear predictors of a logistic fit, projecting onto its principal
# components first for pcaLogistic.
logisticEta <- function(fit, newx) {
    if (!is.null(fit$proj)) {
        newx <- sweep(newx[, fit$proj$genes, drop = FALSE], 2,
            fit$proj$center) %*% fit$proj$rotation
        colnames(newx) <- fit$genes
    }
    sweep(newx[, fit$genes, drop = FALSE] %*% fit$beta, 2, fit$a0, "+")
}

# The predicted class of each cell, as a character vector for every classifier
# (MASS::lda() returns a list with $class, while e1071::svm() and
# randomForest::randomForest() return a factor). `logPrior` holds optional named
# log class priors, added to the classifier's log posterior: the classifier is
# trained on balanced classes, so its own prior is uniform. svm gives no class
# probabilities, so it cannot take priors.
predictClassifierClass <- function(fit, newx, classifier, logPrior = NULL) {
    if (classifier %in% c("logistic", "pcaLogistic")) {
        eta <- logisticEta(fit, newx)
        if (!is.null(logPrior)) {
            eta <- sweep(eta, 2, logPrior[colnames(fit$beta)], "+")
        }
        return(colnames(fit$beta)[max.col(eta, ties.method = "first")])
    }
    if (!is.null(logPrior)) {
        p <- if (classifier == "lda") {
            predict(fit, newx)$posterior
        } else {
            predict(fit, newx, type = "prob")
        }
        p <- sweep(p, 2, exp(logPrior[colnames(p)]), "*")
        return(colnames(p)[max.col(p, ties.method = "first")])
    }
    switch(classifier,
        lda = as.character(predict(fit, newx)$class),
        svm = as.character(predict(fit, newx)),
        rf  = as.character(predict(fit, newx, type = "response"))
    )
}

#' Bootstrap classification confidence for one group
#'
#' Trains a classifier on one group's cells, then repeatedly resamples each
#' cell's counts (multinomial resampling) and classifies the resample. The
#' fraction of resamples assigning a cell to each class is its confidence for
#' that class.
#' [bootstrapLineageConfidence()] runs
#' this on every group of a cell type tree.
#'
#' @param countMat A cells x genes matrix of raw (non-negative integer) counts
#'   for the cells of one group. Normalised or log-transformed values are
#'   rejected, since the resampling needs each cell's total transcript count.
#' @param labels Character vector of the original cell type label for each row
#'   of `countMat`.
#' @param nBoot Number of bootstrap resamples.
#' @param nPerClass Minimum number of training cells per class; classes smaller
#'   than this are up-sampled with replacement.
#' @param cores Number of cores to run the bootstrap resamples on (forked with
#'   BiocParallel); the default runs serially. On Windows, where forking is not
#'   available, it always runs serially. After [set.seed()],
#'   the results are identical whatever the number of cores.
#' @param classifier One of: * `"pcaLogistic"` (the default): ridge
#'   multinomial logistic regression, via glmnet, with the penalty chosen by
#'   cross-validation, on the first 50 principal components of the training
#'   cells (computed exactly for up to 5,000 genes, and with irlba beyond); *
#'   `"logistic"`: the same on all genes; * `"rf"`: random forest, via
#'   randomForest; * `"lda"`: linear discriminant analysis, via MASS; *
#'   `"svm"`: support vector machine, via e1071.
#' @param ntree Number of trees; only used when `classifier = "rf"`.
#' @param priors The classifier is trained with balanced classes, which gives
#'   every class the same prior. `"balanced"` keeps those equal priors;
#'   `"observed"` weights each class by its share of the cells, so a resample
#'   is assigned to its most probable class given how common the classes are.
#'   `"observed"` needs class probabilities, so it is not available for
#'   `classifier = "svm"`.
#' @param correctedCounts `FALSE` (the default) requires whole-number counts.
#'   `TRUE` also accepts counts that are not whole numbers, such as counts
#'   corrected for ambient RNA, and resamples them as they are. They must still
#'   be on the scale of molecule counts: a cell's resample draws as many
#'   transcripts as its total (rounded), so normalised or log-transformed
#'   values would give the wrong amount of resampling noise.
#'
#' @return A numeric matrix with one row per cell (in the order of `countMat`)
#'   and one column per class in `labels`. Entry `[i, k]` is the fraction
#'   of bootstrap resamples in which cell `i` was classified as class `k`,
#'   so each row sums to 1.
#'
#' @examples
#' data(lungSubset)
#' set.seed(1)
#' endo <- lungSubset$cellType %in%
#'     c("Arteriole", "Venous", "Capillary", "Lymphatic")
#' countMat <- t(as.matrix(
#'     SummarizedExperiment::assay(lungSubset, "counts")[, endo]))
#' freq <- bootstrapClassify(countMat, lungSubset$cellType[endo], nBoot = 5)
#' head(freq)
#' @export
bootstrapClassify <- function(countMat, labels, nBoot = 100, nPerClass = 50,
        cores = 1,
        classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"),
        ntree = 500, priors = c("balanced", "observed"),
        correctedCounts = FALSE) {
    classifier <- match.arg(classifier)
    priors <- match.arg(priors)
    checkCores(cores)
    if (!is.logical(correctedCounts) || length(correctedCounts) != 1 ||
        is.na(correctedCounts)) {
        stop("`correctedCounts` must be TRUE or FALSE", call. = FALSE)
    }
    if (classifier == "svm" && priors == "observed") {
        stop("priors = \"observed\" needs class probabilities, which ",
            "classifier = \"svm\" does not give; use another classifier or ",
            "priors = \"balanced\"", call. = FALSE)
    }
    labels <- as.character(labels)
    countMat <- withGeneNames(countMat)
    checkCountMatrix(countMat, labels, correctedCounts)
    checkClassSizes(labels, " to train the classifier")

    blasThreads <- setBlasThreads(1L)
    on.exit(setBlasThreads(blasThreads), add = TRUE)
    features <- classifierFeatures(countMat)
    train <- balancedSample(seq_along(labels), labels, nPerClass)
    fit <- trainClassifier(features$scaled[train, features$keep, drop = FALSE],
        labels[train], classifier, ntree = ntree,
        BPPARAM = parallelParam(cores))

    logPrior <- if (priors == "observed") {
        log(c(table(labels)) / length(labels))
    }
    preds <- bootstrapPredictions(countMat, fit, features$keep, classifier,
        logPrior, nBoot, cores)
    allLabels <- unique(labels)
    freqMatrix <- vapply(allLabels, function(lab) rowMeans(preds == lab),
        numeric(nrow(preds)))
    colnames(freqMatrix) <- allLabels
    freqMatrix
}

# The class predicted for each cell (rows) in each of nBoot resamples of its
# counts (columns). A logistic fit is applied to the resamples' nonzero entries
# alone, without building and scaling each resample as a dense matrix. Each
# resample has its own random stream (see parallelParam()), so the predictions
# do not depend on the number of cores.
bootstrapPredictions <- function(countMat, fit, keep, classifier, logPrior,
        nBoot, cores) {
    layout <- resampleLayout(countMat)
    map <- if (classifier %in% c("logistic", "pcaLogistic")) {
        linearMap(fit, stats::setNames(keep, colnames(countMat)))
    }
    preds <- parLapply(seq_len(nBoot), function(b) {
        if (!is.null(map)) {
            eta <- resampleEta(layout, resampleHits(layout), map)
            if (!is.null(logPrior)) {
                eta <- sweep(eta, 2, logPrior[colnames(eta)], "+")
            }
            return(colnames(eta)[max.col(eta, ties.method = "first")])
        }
        bootMat <- resampleCounts(layout)
        predictClassifierClass(fit,
            scaleCountMat(bootMat)[, keep, drop = FALSE], classifier, logPrior)
    }, parallelParam(cores))
    do.call(cbind, preds)
}

#' Bootstrap confidence scores for every group in a cell type tree
#'
#' For each group in `cellTypeTree` that has children (see
#' [cellTypeTree]), runs [bootstrapClassify()]
#' on that group's cells to get, for every cell, a confidence score for each of
#' the group's child cell types. The result is the `confList` used by the
#' downstream analyses, such as
#' [runProportionsPipeline()] and
#' [runColocalisationModels()]. Cell types
#' without a group (ungrouped leaves) are not refined and keep hard labels.
#'
#' @inheritParams bootstrapClassify
#' @param spe A `SummarizedExperiment` (e.g. `SpatialExperiment`) with
#'   an assay of raw counts, and each cell's type label and unique ID in
#'   `colData`.
#' @param cellType Name of the `colData(spe)` column holding each cell's cell
#'   type label.
#' @param cellID Name of the `colData(spe)` column holding unique cell IDs.
#'   They are returned in the `cellID` column of each result.
#' @param assayName Name of the assay with the raw counts.
#' @param cellTypeTree A cell type tree, see [cellTypeTree].
#' @param skipGroups Names of groups to leave out of the result.
#'
#' @return A list with one data frame per group, named by group. Each has one
#'   row per cell of the group: a `cellID` column, and one column per child cell
#'   type with the cell's confidence in it (between 0 and 1; each row sums to
#'   1). A `"settings"` attribute records how it was computed (`nBoot`,
#'   `nPerClass`, `classifier`, `priors`, `correctedCounts`,
#'   `ntree`, the group's `children` and the `fuzzyCellVersion`).
#'   Groups with no cells, or with only one cell type present, cannot be
#'   bootstrapped: they are skipped with a message.
#'
#' @seealso [crossfitLineageProbabilities()] for
#'   cross-fitted classifier probabilities, an alternative to the bootstrap. The
#'   bootstrap measures how stable each cell's label is when its counts are
#'   resampled; cross-fitting gives each cell the probability of each subtype
#'   from a classifier that did not see the cell, and trains fewer classifiers.
#'   Either result can be passed as `confList` to the downstream analyses.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' set.seed(1)
#' confList <- bootstrapLineageConfidence(lungSubset, lungCellTypeTree,
#'     nBoot = 5)
#' names(confList)
#' head(confList[["Endothelial lineage"]])
#'
#' # the same call with nBoot = 20 gives the confidences in `lungConfList`
#' @export
bootstrapLineageConfidence <- function(spe, cellTypeTree, nBoot = 100,
        nPerClass = 50, cores = 1, skipGroups = character(0),
        classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"),
        ntree = 500, priors = c("balanced", "observed"),
        correctedCounts = FALSE, cellType = "cellType", cellID = "cellID",
        assayName = "counts") {
    classifier <- match.arg(classifier)
    priors <- match.arg(priors)
    checkCores(cores)
    settings <- list(nBoot = nBoot, nPerClass = nPerClass,
        classifier = classifier, priors = priors,
        correctedCounts = correctedCounts, ntree = ntree)
    classify <- function(countMat, labels) {
        bootstrapClassify(countMat, labels, nBoot = nBoot,
            nPerClass = nPerClass, cores = cores, classifier = classifier,
            ntree = ntree, priors = priors, correctedCounts = correctedCounts)
    }
    confidencePerGroup(spe, cellTypeTree, classify, skipGroups, cellType,
        cellID, assayName, method = "bootstrapped", settings = settings)
}

# Runs `classify` on the cells of every group of the tree, except
# `skipGroups`, and returns a list with one data frame per group, named by
# group: a cellID column and one column per child cell type. `classify` takes
# a group's cells x genes count matrix and its labels, and returns a cells x
# child cell types matrix. Used by bootstrapLineageConfidence() and
# crossfitLineageProbabilities(); `method` ("bootstrapped" or
# "cross-fitted") is used in the error raised when no group is left, and
# `settings`, if given, is recorded as an attribute of each table.
confidencePerGroup <- function(spe, cellTypeTree, classify, skipGroups,
        cellType, cellID, assayName, method, settings = NULL) {
    checkSpe(spe)
    checkColumn(spe, cellType, "cellType")
    checkColumn(spe, cellID, "cellID")
    checkAssay(spe, assayName)
    ids <- spe[[cellID]]
    if (anyNA(ids) || anyDuplicated(ids)) {
        example <- if (anyDuplicated(ids)) {
            sprintf(" (e.g. '%s' is repeated)", ids[anyDuplicated(ids)])
        } else {
            ""
        }
        stop("The `cellID` column '", cellID, "' must contain unique, ",
            "non-missing cell IDs", example, call. = FALSE)
    }

    spe <- annotateCellTypeTree(spe, cellTypeTree, cellType = cellType)
    groups <- setdiff(buildCellTypeMap(cellTypeTree)$groupsWithChildren,
        skipGroups)
    confList <- lapply(seq_along(groups), function(gi) {
        groupConfidence(spe, groups, gi, classify, cellType, cellID,
            assayName, settings)
    })
    names(confList) <- groups

    # drop the groups that were skipped
    confList <- confList[!vapply(confList, is.null, logical(1))]
    if (length(confList) == 0) {
        stop("None of the groups could be ", method, ": each had no cells ",
            "or only one cell type present", call. = FALSE)
    }
    confList
}

# The confidences of the cells of group `groups[gi]` from `classify`, as a
# data.frame with a cellID column and one column per child cell type, with
# `settings` (if given) recorded as an attribute. A group needs cells of at
# least 2 cell types; otherwise it is skipped with a message (NULL is
# returned), and its cells keep their hard labels downstream.
groupConfidence <- function(spe, groups, gi, classify, cellType, cellID,
        assayName, settings) {
    g <- groups[gi]
    progress <- sprintf("[%d/%d] Lineage: %s", gi, length(groups), g)
    idx <- which(spe$group == g)
    labels <- spe[[cellType]][idx]
    nClasses <- length(unique(labels))
    if (length(idx) == 0 || nClasses < 2) {
        reason <- if (length(idx) == 0) {
            "no cells in the data"
        } else {
            sprintf("only 1 cell type present ('%s'), needs at least 2",
                labels[1])
        }
        message(progress, " - skipping: ", reason)
        return(NULL)
    }

    message(sprintf("%s (%d cells, %d classes)", progress, length(idx),
        nClasses))
    startTime <- Sys.time()
    confidences <- classify(denseCounts(assay(spe[, idx], assayName)), labels)
    message(sprintf("%s done in %.1fs", progress,
        as.numeric(Sys.time() - startTime, units = "secs")))

    result <- data.frame(cellID = spe[[cellID]][idx], confidences,
        check.names = FALSE)
    if (!is.null(settings)) {
        attr(result, "settings") <- c(settings, list(
            children = unique(labels),
            fuzzyCellVersion = tryCatch(
                as.character(utils::packageVersion("fuzzyCell")),
                error = function(e) NA_character_)))
    }
    result
}

# A cells x genes dense matrix from a genes x cells assay (dense or sparse).
denseCounts <- function(m) {
    if (methods::is(m, "sparseMatrix")) m <- Matrix::as.matrix(m)
    t(as.matrix(m))
}
