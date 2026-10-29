# Cell type tree for the heart transplant data

The cell type hierarchy used with
[`loadAmancherla()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/loadAmancherla.md):
broad lineages, groups of closely related subtypes (for example the T
cell subtypes), and ungrouped leaf cell types. See
[cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md)
for the format.

## Usage

``` r
data(heartCellTypeTree)
```

## Format

A nested list; see
[cellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/cellTypeTree.md).

## Source

Defined in `data-raw/heartCellTypeTree.R`.

## Examples

``` r
data(heartCellTypeTree)
buildCellTypeMap(heartCellTypeTree)$groupsWithChildren
#> [1] "T cells"            "B lineage"          "DC lineage"        
#> [4] "Macrophage lineage" "Fibroblast lineage" "Vascular"          
```
