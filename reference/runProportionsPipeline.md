# Cell type proportion analysis

Computes per-image proportions of every cell type using original labels
(unweighted) and bootstrap confidences (weighted), and tests each
against an outcome with a linear model.
`runProportionsRelativeToParentPipeline()` instead expresses each
subtype as a proportion of its parent group.

## Usage

``` r
runProportionsPipeline(
  spe,
  cellTypeTree,
  confList,
  outcomeCol,
  skipGroups = character(0),
  subject = NULL,
  cellType = "cellType",
  cellID = "cellID",
  imageID = "imageID"
)

runProportionsRelativeToParentPipeline(
  spe,
  cellTypeTree,
  confList,
  outcomeCol,
  skipGroups = character(0),
  subject = NULL,
  cellType = "cellType",
  cellID = "cellID",
  imageID = "imageID"
)
```

## Arguments

- spe:

  A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a cell ID,
  cell type and image for every cell in `colData`.

- cellTypeTree:

  A cell type tree, see
  [cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md).

- confList:

  Confidence list from
  [`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md):
  a named list with one data frame per group, whose names must be groups
  of `cellTypeTree`.

- outcomeCol:

  Name of the per-image outcome column in `colData(spe)`. It must have
  exactly two levels (subset your data to two conditions first) and a
  single value per image: repeated measures and mixed models are not
  supported. An error is raised otherwise, before any computation.

- skipGroups:

  Names of groups of `cellTypeTree` to exclude, as a character vector
  (the default `character(0)` excludes none). An error is raised for a
  name that is not a group of the tree, so the same `skipGroups` can be
  used here and in
  [`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md).

- subject:

  Name of the `colData(spe)` column identifying the subject (for example
  the patient) each image comes from, or `NULL`. If given, every image
  must belong to one subject and every subject must contribute only one
  image, as in
  [`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md).
  With `NULL` this cannot be checked, and a message says so.

- cellType, cellID, imageID:

  Names of the `colData(spe)` columns holding each cell's cell type
  label, a unique cell ID (matching the `cellID` column of `confList`),
  and the image it is in.

## Value

A list with `proportions`, `weighted`, `unweighted` (long `data.frame`s;
only `proportions` and not `weighted`/`unweighted` for the
relative-to-parent version), and `modelWeighted` and `modelUnweighted`
(per cell type model results, where `estimate` is the difference for the
second level of the outcome compared with the first).

## Details

In `runProportionsPipeline()` the unweighted results cover every cell
type and also the whole groups that were bootstrapped. Only the subtypes
of those groups have confidences, so for every other cell type (cell
types without a group, and the group-level proportions) the weighted
proportion is the unweighted one: their labels are not uncertain. The
weighted and unweighted models therefore test the same set of cell
types, and their adjusted p-values are comparable.

The `lineage` column of the results is the group a subtype belongs to
(not the broad lineage). For cell types without a group and for the
group-level proportions it is the cell type or group itself.

## Examples

``` r
data(lungSubset)
data(lungCellTypeTree)
data(lungConfList)
res = runProportionsPipeline(lungSubset, lungCellTypeTree, lungConfList,
    outcomeCol = "disease_status", subject = "subject")
#> Computing unweighted proportions once (shared across all subtypes)...
#> Computing weighted proportions...
#> Fitting weighted models...
head(res$modelWeighted)
#>                      cellType      estimate   statistic      p.value
#> 1           SPP1+ Macrophages -1.310508e-02 -0.77121245 0.4627488770
#> 2    Interstitial Macrophages  3.940842e-02  1.21029437 0.2607089576
#> 3              Monocytes/MDMs -4.222446e-02 -3.99016928 0.0040043674
#> 4 Macrophages - IFN-activated -4.393489e-05 -0.01618578 0.9874825679
#> 5        Alveolar Macrophages -1.883044e-02 -2.55869696 0.0337152591
#> 6                        pDCs  2.360524e-03  5.18432305 0.0008383089
#>         p.adj            lineage
#> 1 0.548469011 Macrophage lineage
#> 2 0.360990284 Macrophage lineage
#> 3 0.021529968 Macrophage lineage
#> 4 0.987482568 Macrophage lineage
#> 5 0.075093077 Macrophage lineage
#> 6 0.008215428    Dendritic cells
data(lungSubset)
data(lungCellTypeTree)
data(lungConfList)
rel = runProportionsRelativeToParentPipeline(lungSubset, lungCellTypeTree,
    lungConfList, outcomeCol = "disease_status", subject = "subject")
#> Computing proportions relative to parent...
#> Fitting unweighted models...
#> Fitting weighted models...
head(rel$modelWeighted)
#>                      cellType     estimate  statistic     p.value       p.adj
#> 1           SPP1+ Macrophages -0.029486426 -0.5204482 0.616844007 0.691612978
#> 2    Interstitial Macrophages  0.385256132  5.0620114 0.000974614 0.006990969
#> 3              Monocytes/MDMs -0.241737210 -3.5653447 0.007343301 0.029427683
#> 4 Macrophages - IFN-activated  0.007970252  0.7049458 0.500849179 0.617713988
#> 5        Alveolar Macrophages -0.122002748 -2.5142228 0.036132460 0.102838540
#> 6                        pDCs  0.160183081  1.3522973 0.213262143 0.343073882
#>              lineage
#> 1 Macrophage lineage
#> 2 Macrophage lineage
#> 3 Macrophage lineage
#> 4 Macrophage lineage
#> 5 Macrophage lineage
#> 6    Dendritic cells
```
