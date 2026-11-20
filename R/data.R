#' A small cropped subset of the pulmonary fibrosis lung data
#'
#' Xenium data of lung biopsies from Vannan et al. (2025), cut down to ship with
#' the package (the full data are `SpatialDatasets::spe_Vannan_2025()`). It
#' has five control and five pulmonary fibrosis biopsies (the five with the
#' fewest cells of each), each cropped to a square window of about 3,000 cells
#' at its centre. Every image comes from a different subject. Cell types with
#' fewer than 20 cells in the windows are removed. The hierarchy for these cell
#' types is in [lungCellTypeTree].
#'
#' @format A `SpatialExperiment` with the raw counts of 343 genes in 29,971
#'   cells (assay `counts`), the cell coordinates in `spatialCoords()`, and
#'   in `colData`:
#' \describe{
#'   \item{cellID}{A unique cell ID.}
#'   \item{imageID}{The image (the window of a biopsy); also `sample_id`.}
#'   \item{subject}{The patient.}
#'   \item{cellType}{The cell type label.}
#'   \item{disease_status}{`"Control"` or `"Disease"`.}
#' }
#' @source Built from `SpatialDatasets::spe_Vannan_2025()` by
#'   `system.file("script", "lungSubset.R", package = "fuzzyCell")`. Original
#'   data: Vannan et al. (2025). Spatial transcriptomics identifies molecular
#'   niche dysregulation associated with distal lung remodeling in pulmonary
#'   fibrosis. Nature Genetics 57, 647-658, \doi{10.1038/s41588-025-02080-x}.
#'   The data are publicly available from NCBI GEO (accession GSE250346); see
#'   the original publication for terms of use.
#' @examples
#' data(lungSubset)
#' lungSubset
#' table(lungSubset$imageID, lungSubset$disease_status)
#' @name lungSubset
#' @docType data
#' @keywords datasets
"lungSubset"

#' Example cell type tree for the lung data
#'
#' An example of the cell type tree a user writes for their own data, here for
#' the cell types of [lungSubset]: broad lineages, groups of closely related
#' subtypes (for example the fibroblast subtypes), and ungrouped cell types.
#' fuzzyCell does not provide trees: for your own data, write a tree like this
#' one for your cell types (see [cellTypeTree] for the format). This tree is
#' included so that the help-page examples can use it.
#'
#' The tree was written by hand by the package authors from the published cell
#' type labels of the Vannan et al. (2025) lung data.
#'
#' @format A nested list; see [cellTypeTree].
#' @source Written in
#'   `system.file("script", "lungCellTypeTree.R", package = "fuzzyCell")`. The
#'   cell type labels are from the original data: Vannan et al. (2025), Nature
#'   Genetics, \doi{10.1038/s41588-025-02080-x}; see [lungSubset].
#' @examples
#' data(lungCellTypeTree)
#' buildCellTypeMap(lungCellTypeTree)$groupsWithChildren
#' @name lungCellTypeTree
#' @docType data
#' @keywords datasets
"lungCellTypeTree"

#' Bootstrap confidences for the lung data
#'
#' The bootstrap classification confidences of the cells in [lungSubset], so
#' that the examples can skip the bootstrap. Made with `set.seed(1)` and then
#' [bootstrapLineageConfidence()] on
#' [lungSubset] and [lungCellTypeTree], with `nBoot = 20` and
#' `classifier = "pcaLogistic"`. Twenty resamples are enough for examples; use
#' at least 100 for real analyses.
#'
#' @format A named list with one data frame per group of [lungCellTypeTree] (see
#'   [bootstrapLineageConfidence()] for the
#'   format).
#' @source Built by
#'   `system.file("script", "lungConfList.R", package = "fuzzyCell")` from
#'   [lungSubset]. Original data: Vannan et al. (2025), Nature Genetics,
#'   \doi{10.1038/s41588-025-02080-x}; see [lungSubset].
#' @examples
#' data(lungConfList)
#' names(lungConfList)
#' head(lungConfList[["Endothelial lineage"]])
#' @name lungConfList
#' @docType data
#' @keywords datasets
"lungConfList"

#' Co-localisation results for the lung data
#'
#' The results of [runColocalisationModels()]
#' for [lungSubset], so that the examples of the plotting functions can skip
#' the analysis. Made with
#' [runColocalisationModels()] on
#' [lungSubset], [lungCellTypeTree] and [lungConfList], with
#' `outcomeCol = "disease_status"` and `subject = "subject"`, and the
#' default radius (50) and window (convex).
#'
#' @format The list returned by
#'   [runColocalisationModels()], without its
#'   `spicy` element (the raw spicyR results), which is left out to keep the
#'   package small.
#' @source Built by
#'   `system.file("script", "lungColocalisation.R", package = "fuzzyCell")` from
#'   [lungSubset], [lungCellTypeTree] and [lungConfList]. Original
#'   data: Vannan et al. (2025), Nature Genetics,
#'   \doi{10.1038/s41588-025-02080-x}; see [lungSubset].
#' @examples
#' data(lungColocalisation)
#' head(lungColocalisation$modelWeighted)
#' @name lungColocalisation
#' @docType data
#' @keywords datasets
"lungColocalisation"
