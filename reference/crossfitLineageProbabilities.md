# Cross-fitted classification probabilities for every group of a tree

An alternative to
[`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md):
for each group with children, a classifier is trained on the group's
cells with balanced classes and gives each cell the probability of each
child cell type. The cells are split into `nFolds` folds and each fold
is scored by a classifier trained on the other folds, so no cell is
scored by a classifier that saw it. Because the classes are balanced in
training, the probabilities give every child cell type the same prior
weight.

## Usage

``` r
crossfitLineageProbabilities(
  spe,
  cellTypeTree,
  nFolds = 5,
  nPerClass = 50,
  seed = 51773,
  classifier = c("pcaLogistic", "logistic", "rf", "lda"),
  ntree = 500,
  priors = c("balanced", "observed"),
  cores = 1,
  skipGroups = character(0),
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

- nFolds:

  Number of folds. The folds are stratified by cell type, so every
  fold's classifier is trained on every child cell type; each child cell
  type therefore needs at least 2 cells.

- nPerClass:

  Minimum number of training cells per class; classes smaller than this
  are up-sampled with replacement.

- seed:

  Random seed, or `NULL` for no fixed seed. Each group is bootstrapped
  with `seed` plus its position in the list of groups, so every group is
  reproducible independently of the others.

- classifier:

  `"pcaLogistic"` (the default), `"logistic"`, `"rf"` or `"lda"`, as in
  [`bootstrapClassify()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapClassify.md).

- ntree:

  Number of trees; only used when `classifier = "rf"`.

- priors:

  `"balanced"` (equal priors, as the classifier is trained) or
  `"observed"` (each class weighted by its share of the group's cells).

- cores:

  Number of cores to fit the folds on in parallel, as in
  [`bootstrapClassify()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapClassify.md).

- skipGroups:

  Names of groups to leave out of the result.

- cellType:

  Name of the `colData(spe)` column holding each cell's cell type label.

- cellID:

  Name of the `colData(spe)` column holding unique cell IDs. They are
  returned in the `cellID` column of each result.

- assayName:

  Name of the assay with the raw counts.

## Value

A confidence list in the format of
[`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md),
with probabilities in place of bootstrap confidences.

## Examples

``` r
data(lungSubset)
data(lungCellTypeTree)
probList = crossfitLineageProbabilities(lungSubset, lungCellTypeTree,
    nFolds = 3, seed = 1)
#> [1/10] Lineage: Macrophage lineage done
#> [2/10] Lineage: Dendritic cells done
#> [3/10] Lineage: Granulocytes done
#> [4/10] Lineage: T cells done
#> [6/10] Lineage: B lineage done
#> [7/10] Lineage: Fibroblast lineage done
#> [8/10] Lineage: Endothelial lineage done
#> [9/10] Lineage: Alveolar done
#> [10/10] Lineage: Airway done
head(probList[["Endothelial lineage"]])
#>               cellID  Arteriole Capillary  Lymphatic     Venous
#> 1 VUHD090_afkejaeg-1 0.09795496 0.5030707 0.37315838 0.02581600
#> 2 VUHD090_afkgbjaj-1 0.38391303 0.3962136 0.11591838 0.10395494
#> 3 VUHD090_afkgcoce-1 0.24468682 0.5278088 0.03473818 0.19276619
#> 4 VUHD090_aflfekgk-1 0.09283809 0.8170683 0.03248653 0.05760710
#> 5 VUHD090_aflfgcca-1 0.04676874 0.8380647 0.04726557 0.06790101
#> 6 VUHD090_aflfkdoo-1 0.16560319 0.3219365 0.18094389 0.33151641
```
