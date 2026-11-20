# Class probabilities from a trained classifier, one row per cell and one
# column per class.
predictClassifierProb <- function(fit, newx, classifier) {
    if (classifier %in% c("logistic", "pcaLogistic")) {
        eta <- logisticEta(fit, newx)
        # subtract each row's largest value before exp() so it cannot overflow
        eta <- eta - eta[cbind(seq_len(nrow(eta)),
            max.col(eta, ties.method = "first"))]
        p <- exp(eta)
        return(p / rowSums(p))
    }
    switch(classifier,
        rf = predict(fit, newx, type = "prob"),
        lda = predict(fit, newx)$posterior
    )
}

# Cross-fitted class probabilities for the cells of one group: the cells are
# split into folds, and each fold is predicted by a classifier trained on the
# other folds with balanced classes (as in bootstrapClassify()), so no cell is
# scored by a classifier that saw it.
crossfitClassify <- function(countMat, labels, nFolds = 5, nPerClass = 50,
        classifier = c("pcaLogistic", "logistic", "rf", "lda"), ntree = 500,
        priors = c("balanced", "observed"), cores = 1) {
    classifier <- match.arg(classifier)
    priors <- match.arg(priors)
    labels <- as.character(labels)
    if (!is.numeric(nFolds) || length(nFolds) != 1 || is.na(nFolds) ||
        nFolds != round(nFolds) || nFolds < 2 || nFolds > length(labels)) {
        stop("`nFolds` must be a whole number between 2 and the number of ",
            "cells (", length(labels), ")", call. = FALSE)
    }

    checkClassSizes(labels, ", so that every fold's training cells include it")
    countMat <- withGeneNames(countMat)
    checkCountMatrix(countMat, labels, correctedCounts = TRUE)
    blasThreads <- setBlasThreads(1L)
    on.exit(setBlasThreads(blasThreads), add = TRUE)
    features <- classifierFeatures(countMat)
    scaledMat <- features$scaled[, features$keep, drop = FALSE]
    classes <- sort(unique(labels))
    fold <- stratifiedFolds(labels, classes, nFolds)
    # one fold per task; each draws its balanced training set from its own
    # random stream
    fits <- parLapply(seq_len(nFolds), function(k) {
        train <- balancedSample(which(fold != k), labels, nPerClass)
        fit <- trainClassifier(scaledMat[train, , drop = FALSE], labels[train],
            classifier, ntree = ntree)
        test <- which(fold == k)
        list(test = test, p = predictClassifierProb(fit,
            scaledMat[test, , drop = FALSE], classifier))
    }, parallelParam(cores))

    probs <- do.call(rbind, lapply(fits, function(f) {
        p <- matrix(0, nrow(f$p), length(classes),
            dimnames = list(NULL, classes))
        p[, colnames(f$p)] <- f$p
        p
    }))
    out <- matrix(0, length(labels), length(classes),
        dimnames = list(NULL, classes))
    out[unlist(lapply(fits, `[[`, "test")), ] <- probs
    if (priors == "observed") {
        shares <- c(table(factor(labels, levels = classes))) / length(labels)
        out <- sweep(out, 2, shares, "*")
        out <- out / rowSums(out)
    }
    out
}

# Each cell's fold: the cells of each class, shuffled, are dealt round the
# folds in turn, so every class with at least 2 cells is in at least 2 folds and
# the fold sizes differ by at most 1.
stratifiedFolds <- function(labels, classes, nFolds) {
    dealt <- unlist(lapply(classes, function(lab) {
        idx <- which(labels == lab)
        idx[sample.int(length(idx))]
    }))
    fold <- integer(length(labels))
    fold[dealt] <- rep(seq_len(nFolds), length.out = length(labels))
    fold
}

#' Cross-fitted classification probabilities for every group of a tree
#'
#' An alternative to
#' [bootstrapLineageConfidence()]: for
#' each group, a classifier is trained on the group's cells with balanced
#' classes and gives each cell the probability of each child cell type. The
#' cells are split into `nFolds` folds and each fold is scored by a classifier
#' trained on the other folds, so no cell is scored by a classifier that saw it.
#' Because the classes are balanced in training, the probabilities give every
#' child cell type the same prior weight.
#'
#' @inheritParams bootstrapLineageConfidence
#' @param nFolds Number of folds. They are stratified by cell type, so each
#'   child cell type needs at least 2 cells.
#' @param classifier `"pcaLogistic"` (the default), `"logistic"`,
#'   `"rf"` or `"lda"`, as in
#'   [bootstrapClassify()].
#' @param priors `"balanced"` (equal priors, as the classifier is trained) or
#'   `"observed"` (each class weighted by its share of the group's cells).
#' @param cores Number of cores to fit the folds on in parallel, as in
#'   [bootstrapClassify()].
#'
#' @return A list with one data frame per group, named by group, in the format
#'   of [bootstrapLineageConfidence()]: a
#'   `cellID` column, and one column per child cell type with the cell's
#'   probability of it (each row sums to 1). Groups with no cells, or with only
#'   one cell type present, are skipped with a message.
#'
#' @seealso [bootstrapLineageConfidence()], the
#'   bootstrap confidences, which measure how stable each cell's label is when
#'   its counts are resampled. Either result can be passed as `confList` to the
#'   downstream analyses.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' set.seed(1)
#' probList <- crossfitLineageProbabilities(lungSubset, lungCellTypeTree,
#'     nFolds = 3)
#' head(probList[["Endothelial lineage"]])
#' @export
crossfitLineageProbabilities <- function(spe, cellTypeTree, nFolds = 5,
        nPerClass = 50, classifier = c("pcaLogistic", "logistic", "rf", "lda"),
        ntree = 500, priors = c("balanced", "observed"), cores = 1,
        skipGroups = character(0), cellType = "cellType", cellID = "cellID",
        assayName = "counts") {
    classifier <- match.arg(classifier)
    priors <- match.arg(priors)
    checkCores(cores)
    classify <- function(countMat, labels) {
        crossfitClassify(countMat, labels, nFolds = nFolds,
            nPerClass = nPerClass, classifier = classifier, ntree = ntree,
            priors = priors, cores = cores)
    }
    confidencePerGroup(spe, cellTypeTree, classify, skipGroups, cellType,
        cellID, assayName, method = "cross-fitted")
}
