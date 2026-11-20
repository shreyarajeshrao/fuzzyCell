## The full-data tests use the real heart transplant data, which
## SpatialDatasets downloads from ExperimentHub. They run only when the
## FUZZYCELL_FULL_TESTS environment variable is set, for example with
## Sys.setenv(FUZZYCELL_FULL_TESTS = "true"), and are skipped otherwise.
skip_unless_full <- function() {
    skip_if(!nzchar(Sys.getenv("FUZZYCELL_FULL_TESTS")),
        "set FUZZYCELL_FULL_TESTS=true to run the full-data tests")
    skip_if_not_installed("SpatialDatasets")
}

## the cell type tree of the heart transplant data (SpatialDatasets::
## spe_Amancherla_2025()), used by the tree tests and the full-data tests
heartCellTypeTree <- list(
    "Immune" = list(
        "T cells" = list(children = c("CD8+ T cells", "CD4+ T cells",
            "Proliferating T cells", "Treg")),
        "B lineage" = list(children = c("B cells", "Plasma")),
        "DC lineage" = list(children = c("cDC1", "cDC2", "pDC",
            "Proliferating DCs", "mDC")),
        "Macrophage lineage" = list(children = c("Macrophages",
            "SPP1+ Macrophages")),
        "NK",
        "Mast"),
    "Mesenchymal" = list(
        "Fibroblast lineage" = list(children = c("Fibroblasts",
            "POSTN+ Fibroblasts", "Myofibroblasts")),
        "Vascular" = list(children = c("Pericytes", "Proliferating pericytes",
            "vSMCs")),
        "Cardiomyocytes",
        "Adipocytes"),
    "Endothelial lineage" = c("Endothelial", "Activated endothelial",
        "BMX+ Activated endothelial", "Proliferating endothelial",
        "Lymphatic endothelial"))

## the pre-transplant heart biopsies of patients with cellular rejection, with
## the columns fuzzyCell uses, and one randomly chosen image per patient
heartPreData <- function() {
    spe <- suppressMessages(SpatialDatasets::spe_Amancherla_2025())
    spe$cellType <- as.character(spe$ct_second_pass)
    spe$x <- spe$adj_x_centroid
    spe$y <- spe$adj_y_centroid
    spe$imageID <- as.character(spe$Sample)
    spe$subject <- as.character(spe$patient_id)
    spe$cellID <- spe$cell_id
    spe <- spe[, !is.na(spe$resolution_broad) &
        spe$patient_rejection_type == "cellular_rejection" &
        spe$biopsy_timing == "pre"]
    images <- withr::with_seed(51773, {
        bySubject <- split(spe$imageID, spe$subject)
        vapply(bySubject, function(im) im[sample.int(length(im), 1L)], "")
    })
    spe[, spe$imageID %in% images]
}

## load the example data and confidence scores once per session
.fullData <- new.env()
fullTestData <- function() {
    if (is.null(.fullData$spe)) {
        .fullData$heartCellTypeTree <- heartCellTypeTree
        .fullData$spe <- heartPreData()
        .fullData$confList <- suppressMessages(suppressWarnings(
            bootstrapLineageConfidence(.fullData$spe,
                .fullData$heartCellTypeTree, nBoot = 5)))
    }
    .fullData
}

## tiny hand-built SpatialExperiment with coordinates in spatialCoords()
miniSpe <- function() {
    counts <- matrix(1:12, nrow = 3,
        dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        colData = data.frame(cellID = paste0("c", 1:4),
            imageID = c("a", "a", "b", "b"),
            cellType = c("T", "T", "B", "B")),
        spatialCoords = cbind(x = c(1, 2, 3, 4), y = c(4, 3, 2, 1)))
}

## small hand-made dataset (6 images, 3 cell types, one group) with
## bootstrap confidences, used to exercise the pipelines end to end
miniPipelineData <- function() {
    set.seed(3)
    nImg <- 6
    nCell <- 70
    cd <- do.call(rbind, lapply(seq_len(nImg), function(i) data.frame(
        cellID = paste0("i", i, "_", seq_len(nCell)),
        imageID = paste0("img", i),
        cellType = sample(c("A", "B", "C"), nCell, TRUE),
        outcome = rep(c("x", "y"), each = 3)[i],
        x = runif(nCell, 0, 300), y = runif(nCell, 0, 300))))
    mu <- matrix(1, 8, nrow(cd))
    mu[1:3, cd$cellType == "A"] <- 5
    mu[4:6, cd$cellType == "B"] <- 5
    counts <- matrix(rpois(length(mu), mu), 8, nrow(cd),
        dimnames = list(paste0("g", 1:8), cd$cellID))
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        colData = cd[, c("cellID", "imageID", "cellType", "outcome")],
        spatialCoords = as.matrix(cd[, c("x", "y")]))
    tree <- list(G = list(children = c("A", "B")), Other = "C")
    conf <- suppressWarnings(suppressMessages(
        bootstrapLineageConfidence(spe, tree, nBoot = 3)))
    list(spe = spe, tree = tree, conf = conf)
}

## a copy of a SpatialExperiment whose cell ID, cell type, image and coordinate
## columns have non-standard names (barcode, ct, sample, px, py)
renamedSpe <- function(spe) {
    spe$barcode <- spe$cellID
    spe$ct <- spe$cellType
    spe$sample <- spe$imageID
    spe$px <- SpatialExperiment::spatialCoords(spe)[, 1]
    spe$py <- SpatialExperiment::spatialCoords(spe)[, 2]
    spe$cellID <- NULL
    spe$cellType <- NULL
    spe$imageID <- NULL
    spe
}

## wrap a per-cell data frame into a SummarizedExperiment (one column per cell),
## for the functions that take an `spe`
asSpe <- function(cells) {
    ids <- paste0("c", seq_len(nrow(cells)))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(0, 1, nrow(cells), dimnames = list("g",
            ids))),
        colData = data.frame(cells, row.names = ids, check.names = FALSE))
}

## hides the messages and warnings of `x`
quiet <- function(x) suppressWarnings(suppressMessages(x))

## a confidence list for group G (cell types A and B) of the cells of `spe`
## with confidence 1 in each cell's own label, so the weighted analyses
## should equal the unweighted ones
hardConf <- function(spe) {
    cells <- getCellData(spe)
    g <- cells[cells$cellType %in% c("A", "B"), ]
    list(G = data.frame(cellID = g$cellID, A = as.numeric(g$cellType == "A"),
        B = as.numeric(g$cellType == "B")))
}
