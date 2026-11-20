## lungConfList: bootstrap confidences for lungSubset, so the examples can skip
## the bootstrap; the same call as in the vignette.
## Source: Vannan et al. (2025). Spatial transcriptomics identifies molecular
## niche dysregulation associated with distal lung remodeling in pulmonary
## fibrosis. Nature Genetics 57, 647-658.
## https://doi.org/10.1038/s41588-025-02080-x
## Data: NCBI GEO accession GSE250346, distributed as a SpatialExperiment by
## SpatialDatasets::spe_Vannan_2025().
## Licence: the data are publicly available from NCBI GEO (GSE250346); see the
## original publication for terms of use.
## Re-create with `source(system.file("script", "lungConfList.R",
## package = "fuzzyCell"))` from the package root, after lungSubset.R and
## lungCellTypeTree.R.

library(fuzzyCell)

load(file.path("data", "lungSubset.rda"))
load(file.path("data", "lungCellTypeTree.rda"))

set.seed(1)
lungConfList <- bootstrapLineageConfidence(lungSubset, lungCellTypeTree,
    nBoot = 20, classifier = "pcaLogistic")

save(lungConfList, file = file.path("data", "lungConfList.rda"),
    compress = "xz", version = 2)
