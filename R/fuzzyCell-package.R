#' fuzzyCell: bootstrap classification confidence for spatial cell types
#'
#' Cell type labels in spatial single-cell data are uncertain, particularly
#' between closely related subtypes. fuzzyCell takes cells that are already
#' labelled and estimates, for each cell, how confident that label is. The
#' confidences are then used as soft cell type weights in downstream analyses,
#' each compared with the same analysis using the original labels.
#'
#' @details
#' A typical analysis has four steps:
#'
#' 1. **Write a cell type tree** for your data: a nested list saying which
#'    subtypes are closely related and may be confused with each other (see
#'    [cellTypeTree]). fuzzyCell does not provide trees; [lungCellTypeTree]
#'    is
#'    an example for the bundled [lungSubset] data. 2. **Estimate
#'    confidences** for every group of related subtypes with
#'    [bootstrapLineageConfidence()] (or, as
#'    an alternative,
#'    [crossfitLineageProbabilities()]). 3.
#'    **Run the analyses**, weighted and unweighted, against a two-level
#'    outcome: cell type proportions
#'    ([runProportionsPipeline()],
#'    [runProportionsRelativeToParentPipeline()]),
#'    co-localisation
#'    ([runColocalisationModels()]) and niches
#'    ([runNicheIdentification()],
#'    [testNicheOutcome()]).
#'    [confidenceWeights()] gives the weights for use
#'    in other tools, such as spicyR and lisaClust. 4. **Plot** the confidences
#'    and results, for example with
#'    [plotWeightedSpatial()] and
#'    [plotProportionBoxplot()].
#'
#' See `vignette("fuzzyCell")` for a full walkthrough.
#'
#' @importFrom stats kmeans lm p.adjust predict as.formula setNames
#' @importFrom rlang .data
#' @importFrom SummarizedExperiment colData assay
#' @import ggplot2
"_PACKAGE"
