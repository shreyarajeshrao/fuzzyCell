# Spatial plots of cell types and confidence scores

Plots the cells of one image, highlighting a cell type or group, or
showing bootstrap confidence scores.

## Usage

``` r
plotWeightedSpatialPair(
  spe,
  cellTypeTree,
  confList,
  imageID,
  groupA,
  subtypeA,
  groupB,
  subtypeB,
  encoding = c("alpha", "colour"),
  cellSize = 1
)

plotHardSpatialPair(
  spe,
  imageID,
  subtypeA,
  subtypeB,
  cellTypeTree = NULL,
  cellSize = 1
)

plotWeightedSpatial(
  spe,
  cellTypeTree,
  confList,
  imageID,
  targetGroup,
  targetSubtype,
  cellSize = 1
)

plotOriginalLabelSpatial(spe, cellTypeTree, imageID, targetGroup, cellSize = 1)
```

## Arguments

- spe:

  A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with a cell ID,
  cell type, image and coordinates for every cell in the form given by
  [`getCellData()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/getCellData.md)
  (`cellID`, `cellType`, `imageID`, `x` and `y`).

- cellTypeTree:

  A cell type tree, see
  [cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md).
  Optional for `plotHardSpatialPair()`, where it is only needed to use a
  group name as a label.

- confList:

  Confidence list from
  [`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md).

- imageID:

  The image to plot. An error is raised if it is not an image in the
  data.

- groupA, groupB, targetGroup:

  Group names. An error is raised for a group that has no cells in the
  data.

- subtypeA, subtypeB, targetSubtype:

  Cell types within the groups, which must have confidence scores in
  `confList` (for `plotHardSpatialPair()`, a cell type or a group name
  in the data).

- encoding:

  How to show the confidence scores: `"alpha"` (the default), as
  transparency, or `"colour"`, as a colour gradient. Either way each of
  the two cell types has its own legend.

- cellSize:

  Size of the highlighted cells, the ones drawn in colour. The grey
  background cells always have size 0.7.

## Value

A `ggplot` object.

## Examples

``` r
data(lungSubset)
data(lungCellTypeTree)
data(lungConfList)
plotHardSpatialPair(lungSubset, "VUHD113", "Mast", "Proliferating FBs",
    lungCellTypeTree)

plotWeightedSpatial(lungSubset, lungCellTypeTree, lungConfList, "VUHD113",
    "Endothelial lineage", "Arteriole")

plotOriginalLabelSpatial(lungSubset, lungCellTypeTree, "VUHD113",
    "Endothelial lineage")

if (requireNamespace("ggnewscale", quietly = TRUE)) {
    plotWeightedSpatialPair(lungSubset, lungCellTypeTree, lungConfList,
        "VUHD113", "Granulocytes", "Mast", "Fibroblast lineage",
        "Proliferating FBs")
}
```
