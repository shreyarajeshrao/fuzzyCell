# A small cropped subset of the pulmonary fibrosis lung data

Xenium data of lung biopsies from Vannan et al. (2025), available in
full as `SpatialDatasets::spe_Vannan_2025()`, cut down to ship with the
package. It has five control and five pulmonary fibrosis biopsies (the
five with the fewest cells of each), each cropped to the roughly 3,000
cells in a square window at its centre, so each image is a window of a
larger biopsy. Every image comes from a different subject. Cell types
with fewer than 20 cells in the windows are removed. The hierarchy for
these cell types is in
[lungCellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/lungCellTypeTree.md).

## Usage

``` r
data(lungSubset)
```

## Format

A `SpatialExperiment` with the raw counts of 343 genes in 29,971 cells
(assay `counts`), the cell coordinates in `spatialCoords()`, and in
`colData`:

- cellID:

  A unique cell ID.

- imageID:

  The image (the window of a biopsy); also `sample_id`.

- subject:

  The patient.

- cellType:

  The cell type label.

- disease_status:

  `"Control"` or `"Disease"`.

## Source

Built from `SpatialDatasets::spe_Vannan_2025()` by
`data-raw/lungSubset.R`.

## Examples

``` r
data(lungSubset)
lungSubset
#> class: SpatialExperiment 
#> dim: 343 29971 
#> metadata(0):
#> assays(1): counts
#> rownames(343): ABCC2 ACKR1 ... YAP1 ZEB1
#> rowData names(0):
#> colnames(29971): TILD167LA_aiflfbgd-1 TILD167LA_aifmhkpg-1 ...
#>   VUILD105MA2_aaaafoki-1 VUILD105MA2_aaaafokj-1
#> colData names(6): cellID imageID ... disease_status sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : x y
#> imgData names(0):
table(lungSubset$imageID, lungSubset$disease_status)
#>              
#>               Control Disease
#>   TILD049MA         0    2989
#>   TILD113LA         0    2997
#>   TILD167LA         0    2999
#>   VUHD038        2997       0
#>   VUHD069        2999       0
#>   VUHD090        2997       0
#>   VUHD095        3000       0
#>   VUHD113        3000       0
#>   VUILD105MA2       0    2996
#>   VUILD142MA        0    2997
```
