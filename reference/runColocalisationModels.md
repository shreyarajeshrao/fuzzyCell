# Co-localisation analysis with soft cell type weights

Tests, for every pair of leaf cell types and groups, whether the
co-localisation of the two differs between the two conditions of an
outcome, once with the original labels and once with the bootstrap
confidences as soft cell type weights. The co-localisation of each image
is the cross-type L-function, computed by spicyR
(`spicyR::spicy(method = "image")`), and relative co-localisation is
Statial's Kontextual test
([`spicyR::kontextualEngine()`](https://sydneybiox.github.io/spicyR/reference/kontextualEngine.html)).

## Usage

``` r
runColocalisationModels(
  spe,
  cellTypeTree,
  confList,
  outcomeCol,
  r = 50,
  subject = NULL,
  cores = 1,
  parentOverrides = character(0),
  autoOverrideParents = TRUE,
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

- outcomeCol:

  Name of the `colData(spe)` column with the per-image outcome. It must
  have exactly two levels and a single value per image.

- r:

  Radius within which cells count as neighbours, in the units of the
  coordinates; a single positive number.

- subject:

  Name of the `colData(spe)` column identifying the subject (for example
  the patient) each image comes from. With `NULL` every image is its own
  subject.

- cores:

  Number of threads.

- parentOverrides:

  Named character vector mapping a `to` label to the ancestor node (a
  group, wrapper or lineage in `cellTypeTree`) that relative
  co-localisation should compare against instead of the label's own
  group, for example `c("Fibroblast lineage" = "Mesenchymal")`.

- autoOverrideParents:

  Whether to add overrides found by
  [`autoParentOverrides()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md).
  Use `FALSE` to condition only on the parents you list in
  `parentOverrides`.

- cellType, cellID, imageID:

  Names of the `colData(spe)` columns holding each cell's cell type
  label, a unique cell ID (matching the `cellID` column of `confList`),
  and the image it is in.

- x, y:

  Names of the columns holding the cell coordinates. For a
  `SpatialExperiment` whose coordinates are in `spatialCoords()`, the
  defaults `"x"` and `"y"` refer to its first two coordinates.

## Value

A list with the tests `modelWeighted`, `modelUnweighted`,
`modelRelativeWeighted` and `modelRelativeUnweighted` (one row per pair:
`estimate`, the difference for the second level of the outcome compared
with the first, its `se`, `statistic`, `p.value` and BH-adjusted
`p.adj`, and for relative co-localisation the `parent` and `side`,
`"attract"` or `"avoid"`), the per-image values `weighted`,
`unweighted`, `relativeWeighted` and `relativeUnweighted` (columns
`imageID`, `from`, `to`, `colocalisation`, and `parent` for the relative
version), and the spicyR results themselves (`spicy`, a list of
`weighted`, `unweighted`, `relativeWeighted` and `relativeUnweighted`),
for spicyR's plotting functions.

## Details

Co-localisation of a pair `from` -\> `to` in an image is the cross-type
L-function: how strongly the `from` cells are attracted to (positive
values) or avoid (negative values) the `to` cells within radius `r`,
with 0 meaning no attraction or repulsion. The area of an image is that
of the square around its cells, and cells near the border are edge
corrected. Each pair is then tested by a linear model of its per-image
values on the outcome (the images are the units); if some subjects have
several images, spicyR fits a mixed model with the subject as a random
effect instead.

With weights (see
[`confidenceWeights()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/confidenceWeights.md)),
a cell of a bootstrapped group belongs to each subtype of the group with
its bootstrap confidence: each cell contributes to a label in proportion
to its weight, and the density of a label is its summed weight per unit
area. With the original labels (weights of 0 and 1) this is the
unweighted L-function.

Relative co-localisation compares the co-localisation of `from` with
`to` to that of `from` with a parent population of `to`, so that
attraction that only arises because both cell types prefer the same
region is not counted. By default the parent of a `to` cell type is its
own group; `parentOverrides` and `autoOverrideParents` give broader
parents. Pairs whose `to` has no parent other than itself are not
tested.

A label is not paired with itself, nor with a group that contains it.

## Examples

``` r
# \donttest{
data(lungSubset)
data(lungCellTypeTree)
data(lungConfList)
res = runColocalisationModels(lungSubset, lungCellTypeTree, lungConfList,
    outcomeCol = "disease_status", subject = "subject")
#> Generating comparisons...
#> Testing weighted co-localisation...
#> Coercing condition into factor. Dropping unused levels. Using disease_status = Control as base comparison group. If this is not the desired base group, please convert cells$disease_status into a factor and change the order of levels(cells$disease_status) so that the base group is at index 1.
#> Warning: Your specified subject parameter has a one-to-one mapping with imageID. Converting to a linear model instead of mixed model.
#> Testing weighted relative co-localisation...
#> Testing unweighted co-localisation...
#> Coercing condition into factor. Dropping unused levels. Using disease_status = Control as base comparison group. If this is not the desired base group, please convert cells$disease_status into a factor and change the order of levels(cells$disease_status) so that the base group is at index 1.
#> Warning: Your specified subject parameter has a one-to-one mapping with imageID. Converting to a linear model instead of mixed model.
#> Warning: essentially perfect fit: summary may be unreliable
#> Warning: essentially perfect fit: summary may be unreliable
#> Warning: essentially perfect fit: summary may be unreliable
#> Warning: essentially perfect fit: summary may be unreliable
#> Warning: essentially perfect fit: summary may be unreliable
#> Warning: essentially perfect fit: summary may be unreliable
#> Testing unweighted relative co-localisation...
head(res$modelWeighted)
#>                  from                 to   estimate       se statistic
#> 1     Dendritic cells Macrophage lineage -22.428813 9.760309 -2.297961
#> 2        Granulocytes Macrophage lineage -11.370974 4.652782 -2.443908
#> 3             T cells Macrophage lineage  -9.451636 7.148071 -1.322264
#> 4           B lineage Macrophage lineage -10.781636 7.197152 -1.498042
#> 5  Fibroblast lineage Macrophage lineage -11.865317 3.981074 -2.980431
#> 6 Endothelial lineage Macrophage lineage  -8.969674 3.893103 -2.303991
#>      p.value     p.adj
#> 1 0.05063170 0.2833686
#> 2 0.04031796 0.2696560
#> 3 0.22263224 0.4833199
#> 4 0.17250037 0.4254510
#> 5 0.01758883 0.2150838
#> 6 0.05015738 0.2833686
# }
```
