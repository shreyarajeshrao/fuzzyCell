# Cell type weights from bootstrap confidences

Builds the soft cell type membership of every cell, the weights that the
weighted analyses pass to spicyR and lisaClust: one column per leaf cell
type and per group of the tree. A cell of a bootstrapped group has, for
each subtype of the group, its bootstrap confidence in that subtype;
every other cell has weight 1 for its own cell type. A group's column is
1 for each of its cells, since the bootstrap reclassifies cells only
within their group.

## Usage

``` r
confidenceWeights(
  spe,
  cellTypeTree,
  confList,
  weighted = TRUE,
  cellType = "cellType",
  cellID = "cellID"
)
```

## Arguments

- spe:

  A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a cell ID,
  cell type, image and coordinates for every cell (see
  [`getCellData()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/getCellData.md)).

- cellTypeTree:

  A cell type tree, see
  [cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md).

- confList:

  Confidence list from
  [`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md).

- weighted:

  `TRUE` for the bootstrap confidences, `FALSE` for the original labels
  (0/1 weights) on the same columns.

- cellType, cellID:

  Names of the `colData(spe)` columns holding each cell's cell type
  label and a unique cell ID (matching the `cellID` column of
  `confList`).

## Value

A numeric matrix with one row per cell of `spe` (in its order) and one
named column per leaf cell type and group.

## Examples

``` r
data(lungSubset)
data(lungCellTypeTree)
data(lungConfList)
w = confidenceWeights(lungSubset, lungCellTypeTree, lungConfList)
head(w[, c("Arteriole", "Venous", "Endothelial lineage")])
#>      Arteriole Venous Endothelial lineage
#> [1,]         0      0                   0
#> [2,]         0      0                   0
#> [3,]         0      0                   0
#> [4,]         0      0                   0
#> [5,]         0      0                   0
#> [6,]         0      0                   0
```
