# Cell type tree for the lung data

The cell type hierarchy for
[lungSubset](https://shreyarajeshrao.github.io/fuzzyCell/reference/lungSubset.md)
(and the full `SpatialDatasets::spe_Vannan_2025()` data): broad
lineages, groups of closely related subtypes (for example the fibroblast
subtypes), and ungrouped leaf cell types. See
[cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md)
for the format.

## Usage

``` r
data(lungCellTypeTree)
```

## Format

A nested list; see
[cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md).

## Source

Defined in `data-raw/lungCellTypeTree.R`.

## Examples

``` r
data(lungCellTypeTree)
buildCellTypeMap(lungCellTypeTree)$groupsWithChildren
#>  [1] "Macrophage lineage"  "Dendritic cells"     "Granulocytes"       
#>  [4] "T cells"             "NK lineage"          "B lineage"          
#>  [7] "Fibroblast lineage"  "Endothelial lineage" "Alveolar"           
#> [10] "Airway"             
```
