# Bootstrap confidence scores for every group in a cell type tree

For each group in `cellTypeTree` that has children (see
[cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md)),
runs
[`bootstrapClassify()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapClassify.md)
on that group's cells to get, for every cell, a confidence score for
each of the group's child cell types. The result is the `confList` used
by the downstream analyses, such as
[`runProportionsPipeline()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runProportionsPipeline.md)
and
[`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md).
Cell types without a group (ungrouped leaves) are not refined and keep
hard labels.

## Usage

``` r
bootstrapLineageConfidence(
  spe,
  cellTypeTree,
  nBoot = 100,
  nPerClass = 50,
  cores = 1,
  seed = 51773,
  skipGroups = character(0),
  subprocess = FALSE,
  classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"),
  ntree = 500,
  priors = c("balanced", "observed"),
  correctedCounts = FALSE,
  cellType = "cellType",
  cellID = "cellID",
  assayName = "counts"
)
```

## Arguments

- spe:

  A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a counts
  assay and per-cell cell type labels and unique cell IDs in `colData`.

- cellTypeTree:

  A cell type tree, see
  [cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md).

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

  Random seed, or `NULL` for no fixed seed. Each group is bootstrapped
  with `seed` plus its position in the list of groups, so every group is
  reproducible independently of the others.

- skipGroups:

  Names of groups to leave out of the result.

- subprocess:

  If `TRUE`, each group is bootstrapped in a fresh R session (via
  callr), which only holds that group's counts. With `cores` above 1,
  every worker starts as a copy of the session it is forked from, so
  this keeps the workers small on large data. It requires the same
  version of fuzzyCell to be installed.

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

- cellType:

  Name of the `colData(spe)` column holding each cell's cell type label.

- cellID:

  Name of the `colData(spe)` column holding unique cell IDs. They are
  returned in the `cellID` column of each result.

- assayName:

  Name of the assay with the raw counts.

## Value

A named list with one element per group (named by group). Groups with no
cells, or with only one cell type present, cannot be bootstrapped: they
are skipped with a message and left out of the list. Each element is a
`data.frame` with one row per cell in that group, a `cellID` column, and
one column per child cell type of the group holding the cell's
confidence score for that cell type (between 0 and 1, rows summing to 1
across the child columns). Each data frame also carries a `"settings"`
attribute recording how it was computed (`nBoot`, `nPerClass`,
`classifier`, `priors`, `correctedCounts`, `ntree`, the group's `seed`,
its `children` and the `fuzzyCellVersion`).

## Examples

``` r
data(lungSubset)
data(lungCellTypeTree)
confList = bootstrapLineageConfidence(lungSubset, lungCellTypeTree, nBoot = 5,
    seed = 1)
#> [1/10] Lineage: Macrophage lineage (4571 cells, 5 classes)
#> Minimum class size: 50
#> Finished training classifier.
#> [1/10] Lineage: Macrophage lineage done in 0.5s
#> [2/10] Lineage: Dendritic cells (193 cells, 3 classes)
#> Minimum class size: 50
#> Finished training classifier.
#> [2/10] Lineage: Dendritic cells done in 0.1s
#> [3/10] Lineage: Granulocytes (1575 cells, 3 classes)
#> Minimum class size: 50
#> Finished training classifier.
#> [3/10] Lineage: Granulocytes done in 0.1s
#> [4/10] Lineage: T cells (2927 cells, 3 classes)
#> Minimum class size: 165
#> Finished training classifier.
#> [4/10] Lineage: T cells done in 0.3s
#> [5/10] Lineage: NK lineage - skipping: only 1 cell type present ('NK/NKT'), needs at least 2
#> [6/10] Lineage: B lineage (1210 cells, 2 classes)
#> Minimum class size: 359
#> Finished training classifier.
#> [6/10] Lineage: B lineage done in 0.3s
#> [7/10] Lineage: Fibroblast lineage (5148 cells, 7 classes)
#> Minimum class size: 50
#> Finished training classifier.
#> [7/10] Lineage: Fibroblast lineage done in 0.5s
#> [8/10] Lineage: Endothelial lineage (6214 cells, 4 classes)
#> Minimum class size: 108
#> Finished training classifier.
#> [8/10] Lineage: Endothelial lineage done in 0.3s
#> [9/10] Lineage: Alveolar (3438 cells, 4 classes)
#> Minimum class size: 50
#> Finished training classifier.
#> [9/10] Lineage: Alveolar done in 0.4s
#> [10/10] Lineage: Airway (1899 cells, 6 classes)
#> Minimum class size: 100
#> Finished training classifier.
#> [10/10] Lineage: Airway done in 0.3s
names(confList)
#> [1] "Macrophage lineage"  "Dendritic cells"     "Granulocytes"       
#> [4] "T cells"             "B lineage"           "Fibroblast lineage" 
#> [7] "Endothelial lineage" "Alveolar"            "Airway"             
head(confList[["Endothelial lineage"]])
#>               cellID Capillary Venous Lymphatic Arteriole
#> 1 VUHD090_afkejaeg-1       0.6    0.0       0.4       0.0
#> 2 VUHD090_afkgbjaj-1       1.0    0.0       0.0       0.0
#> 3 VUHD090_afkgcoce-1       0.8    0.2       0.0       0.0
#> 4 VUHD090_aflfekgk-1       0.8    0.0       0.0       0.2
#> 5 VUHD090_aflfgcca-1       1.0    0.0       0.0       0.0
#> 6 VUHD090_aflfkdoo-1       0.6    0.0       0.2       0.2

# the same call with nBoot = 20 gives the confidences in `lungConfList`
```
