# Bootstrap confidences for the lung data

The bootstrap classification confidences of the cells in
[lungSubset](https://shreyarajeshrao.github.io/fuzzyCell/reference/lungSubset.md),
from
`bootstrapLineageConfidence(lungSubset, lungCellTypeTree, nBoot = 20, classifier = "pcaLogistic", seed = 1)`,
so that the examples can skip the bootstrap.

## Usage

``` r
data(lungConfList)
```

## Format

A named list with one data frame per group of
[lungCellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/lungCellTypeTree.md)
(see
[`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md)
for the format).

## Source

Built by `data-raw/lungConfList.R`.

## Examples

``` r
data(lungConfList)
names(lungConfList)
#> [1] "Macrophage lineage"  "Dendritic cells"     "Granulocytes"       
#> [4] "T cells"             "B lineage"           "Fibroblast lineage" 
#> [7] "Endothelial lineage" "Alveolar"            "Airway"             
head(lungConfList[["Endothelial lineage"]])
#>               cellID Capillary Venous Lymphatic Arteriole
#> 1 VUHD090_afkejaeg-1      0.60   0.00      0.40      0.00
#> 2 VUHD090_afkgbjaj-1      1.00   0.00      0.00      0.00
#> 3 VUHD090_afkgcoce-1      0.90   0.10      0.00      0.00
#> 4 VUHD090_aflfekgk-1      0.95   0.00      0.00      0.05
#> 5 VUHD090_aflfgcca-1      0.85   0.15      0.00      0.00
#> 6 VUHD090_aflfkdoo-1      0.70   0.00      0.05      0.25
```
