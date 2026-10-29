# Heart transplant example data

Downloads the Amancherla et al. (2025) heart transplant Xenium dataset
from the SpatialDatasets package (via ExperimentHub) and prepares it for
fuzzyCell: renames the columns fuzzyCell expects, keeps patients with
cellular rejection, and keeps one randomly chosen image per patient for
each biopsy timing. The first call downloads the data, which is large.
With `"both"` every patient has an image before and after transplant, so
those images are not independent and cannot be passed together to
[`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md)
with `subject = "subject"`.

## Usage

``` r
loadAmancherla(timing = c("pre", "post", "both"), seed = 51773)
```

## Arguments

- timing:

  Which biopsies to keep: `"pre"` or `"post"` transplant, or `"both"`.

- seed:

  Random seed used to choose one image per patient.

## Value

A `SpatialExperiment`.

## Details

The prepared object has the columns fuzzyCell functions use: `cellID`,
`imageID`, `cellType`, `x`, `y` and `subject`, and the per-patient
outcome `resolution_broad` (`"persistence"` or `"resolved"`). The cell
type hierarchy for these cell types is in
[heartCellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/heartCellTypeTree.md).

## Examples

``` r
# \donttest{
spe = loadAmancherla("pre")
#> Warning: previous export ''Ferguson_Images'' is being replaced
#> Warning: previous export ''spe_Ferguson_2022'' is being replaced
#> Warning: previous export ''fergusonClinical'' is being replaced
#> Warning: previous export ''spe_Schurch_2020'' is being replaced
#> Warning: previous export ''spe_Ali_2020'' is being replaced
#> see ?SpatialDatasets and browseVignettes('SpatialDatasets') for documentation
#> loading from cache
spe
#> class: SpatialExperiment 
#> dim: 477 30878 
#> metadata(0):
#> assays(3): counts logcounts scaledata
#> rownames(477): ABCC11 ABHD16A ... XCR1 ZCCHC12
#> rowData names(0):
#> colnames(30878): aaliapfm-1_1 aaljojfh-1_1 ... oimehgjn-1_2
#>   oimeocni-1_2
#> colData names(61): X orig.ident ... subject cellID
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : x y
#> imgData names(1): sample_id
table(spe$resolution_broad)
#> 
#> persistence    resolved 
#>       16676       14202 
# }
```
