## lungCellTypeTree: the cell type hierarchy for lungSubset.
## Written by hand by the package authors from the published cell type labels
## (column final_CT) of the Vannan et al. (2025) pulmonary fibrosis lung Xenium
## data, as an example of a user-supplied tree (see ?cellTypeTree). It covers
## only the cell types present in lungSubset.
## Source: Vannan et al. (2025). Spatial transcriptomics identifies molecular
## niche dysregulation associated with distal lung remodeling in pulmonary
## fibrosis. Nature Genetics 57, 647-658.
## https://doi.org/10.1038/s41588-025-02080-x
## Data: NCBI GEO accession GSE250346, distributed as a SpatialExperiment by
## SpatialDatasets::spe_Vannan_2025().
## Licence: the data are publicly available from NCBI GEO (GSE250346); see the
## original publication for terms of use.
## Re-create with `source(system.file("script", "lungCellTypeTree.R",
## package = "fuzzyCell"))` from the package root.

lungCellTypeTree <- list(
    "Immune" = list(
        "Myeloid" = list(
            "Macrophage lineage" = list(children = c("Alveolar Macrophages",
                "Interstitial Macrophages", "SPP1+ Macrophages",
                "Macrophages - IFN-activated", "Monocytes/MDMs")),
            "Dendritic cells" = list(children = c("cDCs", "Migratory DCs",
                "pDCs")),
            "Granulocytes" = list(children = c("Mast", "Basophils",
                "Neutrophils")),
            "Proliferating Myeloid"),
        "Lymphoid" = list(
            "T cells" = list(children = c("CD4+ T-cells", "CD8+ T-cells",
                "Tregs")),
            "B lineage" = list(children = c("B cells", "Plasma")),
            "NK/NKT")),
    "Mesenchymal" = list(
        "Fibroblast lineage" = list(children = c("Alveolar FBs",
            "Adventitial FBs", "Activated Fibrotic FBs", "Inflammatory FBs",
            "Subpleural FBs", "Proliferating FBs", "Myofibroblasts")),
        "SMCs/Pericytes"),
    "Endothelial lineage" = list(children = c("Arteriole", "Venous",
        "Capillary", "Lymphatic")),
    "Epithelial lineage" = list(
        "Alveolar" = list(children = c("AT1", "AT2", "Transitional AT2",
            "Proliferating AT2")),
        "Airway" = list(children = c("Basal", "Secretory", "Goblet",
            "Multiciliated", "RASC", "KRT5-/KRT17+"))))

save(lungCellTypeTree, file = file.path("data", "lungCellTypeTree.rda"),
    compress = "xz", version = 2)
