# Bootstrap classification confidence for one group

The single-group worker behind bootstrapLineageConfidence(). Trains a
classifier on one group's cells, then repeatedly resamples each cell's
counts (multinomial resampling) and re-predicts its class. The fraction
of resamples assigning a cell to each class is its confidence for that
class.

## Usage

``` r
bootstrapClassify(
  countMat,
  labels,
  nBoot = 100,
  nPerClass = 50,
  cores = 1,
  seed = NULL,
  classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"),
  ntree = 500,
  priors = c("balanced", "observed"),
  correctedCounts = FALSE
)
```

## Arguments

- countMat:

  A cells x genes matrix of raw (non-negative integer) counts for the
  cells of one group. Normalised or log-transformed values are rejected,
  since the resampling needs each cell's total transcript count.

- labels:

  Character vector of the original cell type label for each row of
  `countMat`.

- nBoot:

  Number of bootstrap resamples.

- nPerClass:

  Minimum number of training cells per class; classes smaller than this
  are up-sampled with replacement.

- cores:

  Number of cores to run the bootstrap resamples on (forked with
  BiocParallel); the default runs serially. On Windows, where forking is
  not available, it always runs serially. If `seed` is given, the
  results are identical whatever the number of cores.

- seed:

  Random seed for reproducibility, or `NULL`. Recommended for
  reproducibility.

- classifier:

  One of `"pcaLogistic"` (ridge multinomial logistic regression, via
  glmnet, with the penalty chosen by cross-validation, on the first 50
  principal components of the training cells, computed exactly for up to
  5,000 genes and with irlba beyond; the default), `"logistic"` (the
  same on all genes), `"rf"` (random forest, via randomForest), `"lda"`
  (linear discriminant analysis, via MASS) or `"svm"` (support vector
  machine, via e1071).

- ntree:

  Number of trees; only used when `classifier = "rf"`.

- priors:

  The classifier is trained with balanced classes, which gives every
  class the same prior. `"balanced"` keeps those equal priors;
  `"observed"` weights each class by its share of the cells, so a
  resample is assigned to its most probable class given how common the
  classes are.

- correctedCounts:

  `FALSE` (the default) requires whole-number counts. `TRUE` also
  accepts counts that are not whole numbers, such as counts corrected
  for ambient RNA, and resamples them as they are. They must still be on
  the scale of molecule counts: a cell's resample draws as many
  transcripts as its total (rounded), so normalised or log-transformed
  values would give the wrong amount of resampling noise.

## Value

A numeric matrix with one row per cell (in the order of `countMat`) and
one column per class in `labels`. Entry `[i, k]` is the fraction of
bootstrap resamples in which cell `i` was classified as class `k`, so
each row sums to 1.

## Examples

``` r
data(lungSubset)
endo = lungSubset$cellType %in% c("Arteriole", "Venous", "Capillary", "Lymphatic")
#> Loading required namespace: SpatialExperiment
countMat = t(as.matrix(SummarizedExperiment::assay(lungSubset, "counts")[, endo]))
freq = bootstrapClassify(countMat, lungSubset$cellType[endo], nBoot = 5, seed = 1)
#> Minimum class size: 108
#> Finished training classifier.
head(freq)
#>      Capillary Venous Lymphatic Arteriole
#> [1,]       0.4    0.0       0.6         0
#> [2,]       1.0    0.0       0.0         0
#> [3,]       1.0    0.0       0.0         0
#> [4,]       1.0    0.0       0.0         0
#> [5,]       1.0    0.0       0.0         0
#> [6,]       0.8    0.2       0.0         0
```
