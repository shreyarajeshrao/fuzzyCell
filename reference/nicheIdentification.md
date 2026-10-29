# Niche identification

Defines niches by k-means clustering of per-cell feature vectors, once
with hard cell type membership and once with bootstrap-confidence
weights. The features of each cell are its local indicators of spatial
association from lisaClust
([`lisaClust::lisa()`](https://rdrr.io/pkg/lisaClust/man/lisa.html)) at
several radii.

## Usage

``` r
runNicheIdentification(
  spe,
  cellTypeTree,
  confList,
  toLabels = NULL,
  r = c(20, 50, 100),
  k,
  window = "convex",
  lisaFunc = "K",
  cores = 1,
  seed = 51773,
  cellType = "cellType",
  cellID = "cellID",
  imageID = "imageID",
  x = "x",
  y = "y"
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

- toLabels:

  Labels (leaves and/or groups) to build features from; defaults to
  every leaf and group in the tree. An error is raised for a label that
  is not in the tree.

- r:

  Radii of the local curves: one or more positive numbers.

- k:

  Number of niches, a whole number of at least 2 and below the number of
  cells.

- window, lisaFunc:

  The image window and the local function (`"K"` or `"L"`), as in
  [`lisaClust::lisa()`](https://rdrr.io/pkg/lisaClust/man/lisa.html).

- cores:

  Number of threads.

- seed:

  Random seed for the k-means clustering, as in
  [`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md),
  or `NULL` for no fixed seed. The same seed is used for the unweighted
  and weighted clusterings.

- cellType, cellID, imageID:

  Names of the `colData(spe)` columns holding each cell's cell type
  label, a unique cell ID (matching the `cellID` column of `confList`),
  and the image it is in.

- x, y:

  Names of the columns holding the cell coordinates. For a
  `SpatialExperiment` whose coordinates are in `spatialCoords()`, the
  defaults `"x"` and `"y"` refer to its first two coordinates.

## Value

A list with the feature matrices (`lisaMatUnweighted`,
`lisaMatWeighted`) and named cell-to-niche vectors `nicheUnweighted` and
`nicheWeighted`.

## Details

For every cell, every label in `toLabels` and every radius in `r`, the
feature is the local co-localisation of the cell's neighbours with that
label (lisaClust's local K or L function, centred so that 0 means no
attraction or repulsion). The weighted features count each neighbour
with its bootstrap confidence of having the label, instead of 0 or 1
(see
[`confidenceWeights()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/confidenceWeights.md)),
and a label's density is its summed weight per unit area.

Niche names (`niche1`, `niche2`, ...) are arbitrary. K-means starts from
random cluster centres, so the niches depend on the random seed. The
unweighted and weighted clusterings are both started from the same
`seed`, so that if the weights made no difference the two sets of niches
would be identical.

## Examples

``` r
data(lungSubset)
data(lungCellTypeTree)
data(lungConfList)
niches = runNicheIdentification(lungSubset, lungCellTypeTree, lungConfList,
    r = 50, k = 6)
#> Computing unweighted LISA curves...
#> Generating local indicators of spatial association.
#> Computing weighted LISA curves...
#> Generating local indicators of spatial association.
#> Clustering niches...
table(weighted = niches$nicheWeighted, unweighted = niches$nicheUnweighted)
#>         unweighted
#> weighted niche1 niche2 niche3 niche4 niche5 niche6
#>   niche1  13204     13      0    157      0    164
#>   niche2     19   1946      0     17      1      3
#>   niche3      0      0    396      0      5      0
#>   niche4    194     14      0   8326      0    436
#>   niche5      0      0    166     16    773    426
#>   niche6    233      0      0    109     24   3329
```
