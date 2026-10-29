# Per-cell data frame from a SummarizedExperiment

Returns `colData(spe)` as a `data.frame`. For a `SpatialExperiment`
whose coordinates are stored in `spatialCoords()` rather than in
`colData()`, the coordinates are added as columns `x` and `y`. The
spatial functions in fuzzyCell expect `cellID`, `imageID`, `cellType`,
`x` and `y` columns.

## Usage

``` r
getCellData(spe)
```

## Arguments

- spe:

  A `SummarizedExperiment` or `SpatialExperiment`.

## Value

A `data.frame` with one row per cell.

## Examples

``` r
data(lungSubset)
head(getCellData(lungSubset))
#>                                    cellID   imageID subject  cellType
#> TILD167LA_aiflfbgd-1 TILD167LA_aiflfbgd-1 TILD167LA TILD167     Basal
#> TILD167LA_aifmhkpg-1 TILD167LA_aifmhkpg-1 TILD167LA TILD167      RASC
#> TILD167LA_aifminbe-1 TILD167LA_aifminbe-1 TILD167LA TILD167      RASC
#> TILD167LA_aifogebk-1 TILD167LA_aifogebk-1 TILD167LA TILD167     Basal
#> TILD167LA_aifoiaih-1 TILD167LA_aifoiaih-1 TILD167LA TILD167     Basal
#> TILD167LA_aifpljci-1 TILD167LA_aifpljci-1 TILD167LA TILD167 Secretory
#>                      disease_status sample_id        x         y
#> TILD167LA_aiflfbgd-1        Disease TILD167LA 23675.95 -1982.301
#> TILD167LA_aifmhkpg-1        Disease TILD167LA 23651.80 -2217.509
#> TILD167LA_aifminbe-1        Disease TILD167LA 23666.44 -2218.133
#> TILD167LA_aifogebk-1        Disease TILD167LA 23616.25 -2217.349
#> TILD167LA_aifoiaih-1        Disease TILD167LA 23687.30 -2127.461
#> TILD167LA_aifpljci-1        Disease TILD167LA 23675.64 -1988.052
```
