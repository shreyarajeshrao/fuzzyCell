# Package index

## Cell-type trees

- [`buildCellTypeMap()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md)
  [`annotateCellTypeTree()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md)
  [`resolveTreeLeaves()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md)
  [`autoParentOverrides()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md)
  : Cell type tree utilities

## Classification confidence

- [`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md)
  : Bootstrap confidence scores for every group in a cell type tree
- [`bootstrapClassify()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapClassify.md)
  : Bootstrap classification confidence for one group
- [`crossfitLineageProbabilities()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/crossfitLineageProbabilities.md)
  : Cross-fitted classification probabilities for every group of a tree
- [`confidenceWeights()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/confidenceWeights.md)
  : Cell type weights from bootstrap confidences

## Proportions and co-localisation

- [`runProportionsPipeline()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runProportionsPipeline.md)
  [`runProportionsRelativeToParentPipeline()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runProportionsPipeline.md)
  : Cell type proportion analysis
- [`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md)
  : Co-localisation analysis with soft cell type weights
- [`plotProportionBoxplot()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotProportionBoxplot.md)
  : Plot one cell type's proportion by outcome
- [`plotMetricBoxplot()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotMetricBoxplot.md)
  : Plot one pair's co-localisation by outcome

## Niches

- [`runNicheIdentification()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nicheIdentification.md)
  : Niche identification
- [`testNicheOutcome()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/testNicheOutcome.md)
  : Test niche proportions against an outcome
- [`plotNicheSpatial()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nichePlots.md)
  [`nicheCompositionMap()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nichePlots.md)
  [`nicheBarplot()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nichePlots.md)
  [`nicheProportionBoxplot()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nichePlots.md)
  : Niche plots

## Spatial plots

- [`plotWeightedSpatialPair()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotSpatial.md)
  [`plotHardSpatialPair()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotSpatial.md)
  [`plotWeightedSpatial()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotSpatial.md)
  [`plotOriginalLabelSpatial()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotSpatial.md)
  : Spatial plots of cell types and confidence scores
- [`plotOneMarkerSpatial()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotMarkerSpatial.md)
  [`plotTwoMarkerSpatial()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/plotMarkerSpatial.md)
  : Spatial plots of marker expression within a group

## Data

- [`getCellData()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/getCellData.md)
  : Per-cell data frame from a SummarizedExperiment
- [`loadAmancherla()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/loadAmancherla.md)
  : Heart transplant example data
- [`lungSubset`](https://shreyarajeshrao.github.io/fuzzyCell/reference/lungSubset.md)
  : A small cropped subset of the pulmonary fibrosis lung data
- [`lungCellTypeTree`](https://shreyarajeshrao.github.io/fuzzyCell/reference/lungCellTypeTree.md)
  : Cell type tree for the lung data
- [`lungConfList`](https://shreyarajeshrao.github.io/fuzzyCell/reference/lungConfList.md)
  : Bootstrap confidences for the lung data
- [`heartCellTypeTree`](https://shreyarajeshrao.github.io/fuzzyCell/reference/heartCellTypeTree.md)
  : Cell type tree for the heart transplant data
