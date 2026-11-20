## lungColocalisation: the co-localisation results for lungSubset, so the
## examples of the plotting functions can skip runColocalisationModels(). The
## raw spicyR results (the `spicy` element) are left out to keep the package
## small.
## Source: Vannan et al. (2025). Spatial transcriptomics identifies molecular
## niche dysregulation associated with distal lung remodeling in pulmonary
## fibrosis. Nature Genetics 57, 647-658.
## https://doi.org/10.1038/s41588-025-02080-x
## Data: NCBI GEO accession GSE250346, distributed as a SpatialExperiment by
## SpatialDatasets::spe_Vannan_2025().
## Licence: the data are publicly available from NCBI GEO (GSE250346); see the
## original publication for terms of use.
## Re-create with `source(system.file("script", "lungColocalisation.R",
## package = "fuzzyCell"))` from the package root, after lungSubset.R,
## lungCellTypeTree.R and lungConfList.R.

library(fuzzyCell)

load(file.path("data", "lungSubset.rda"))
load(file.path("data", "lungCellTypeTree.rda"))
load(file.path("data", "lungConfList.rda"))

lungColocalisation <- runColocalisationModels(lungSubset, lungCellTypeTree,
    lungConfList, outcomeCol = "disease_status", subject = "subject")
lungColocalisation$spicy <- NULL

save(lungColocalisation, file = file.path("data", "lungColocalisation.rda"),
    compress = "xz", version = 2)
