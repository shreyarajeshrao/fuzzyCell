# Cell type tree utilities

Helpers for the cell type tree: a nested list describing how cell types
are organised into broad lineages and (optionally) groups of closely
related subtypes.
[`bootstrapLineageConfidence()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/bootstrapLineageConfidence.md)
refines the cells within each group, and the analysis functions compare
cell types, groups and lineages.

## Usage

``` r
buildCellTypeMap(cellTypeTree)

annotateCellTypeTree(spe, cellTypeTree, cellType = "cellType")

resolveTreeLeaves(cellTypeTree, nodeName)

autoParentOverrides(cellTypeTree)
```

## Arguments

- cellTypeTree:

  A nested list describing the cell type hierarchy (see Details and
  [heartCellTypeTree](https://shreyarajeshrao.github.io/fuzzyCell/reference/heartCellTypeTree.md)).

- spe:

  A `SummarizedExperiment` (e.g. a `SpatialExperiment`), or a
  `data.frame`/`DataFrame` of per-cell data.

- cellType:

  Name of the column of `colData(spe)` (or of the data frame) holding
  each cell's cell type label. An error is raised if it is not a column.

- nodeName:

  A single string: the name of any node in `cellTypeTree`, a lineage, an
  organisational wrapper, a group or a cell type.

## Value

- `buildCellTypeMap()`: a list with `cellMap` (a `data.frame` with
  columns `cellType`, `group` and `broad`, one row per cell type) and
  `groupsWithChildren` (names of the groups to be bootstrap-refined). An
  error is raised if the tree is not a valid tree as described in
  Details: it must be a non-empty named list, every top-level entry must
  be named, groups may only contain `children`, cell type names must be
  non-empty strings, and names must not be repeated or shared between
  cell types and lineages, wrappers or groups.

- `annotateCellTypeTree()`: `spe` with `group` and `lineage` columns
  added (a message reports any existing `group` or `lineage` column
  whose values are replaced). An error is raised if any cell's cell type
  is not in the tree (including `NA`).

- `resolveTreeLeaves()`: character vector of every cell type beneath
  `nodeName`.

- `autoParentOverrides()`: named character vector mapping each group
  that has a named ancestor to that ancestor, for use as
  `parentOverrides` in
  [`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md).

## Details

A tree is a named list with one entry per broad lineage (these become
the `lineage` column). Each lineage is either a character vector of cell
types, or a list whose entries are:

- a bare string: a cell type that is not in any group;

- a named `list(children = c(...))`: a **group** of closely related cell
  types, the unit that is bootstrap-refined. A group may contain only
  `children`;

- a named list of further entries: purely organisational nesting (a
  "wrapper") that is not a group itself, to any depth.

A lineage can also itself be `list(children = c(...))`, in which case
the lineage name doubles as the group name. A cell type's group is its
nearest `list(children = ...)` ancestor, however many wrappers sit above
it. In the example below `CD8_Resting` is in the group `CD8+ T cells`,
not `T cells`.

Every cell type must appear exactly once, and lineage, wrapper and group
names must be unique and different from every cell type name.

    list(
        "Immune" = list(
            "T cells" = list(
                "CD8+ T cells" = list(children = c("CD8_Resting", "CD8_Memory")),
                "CD4+ T cells" = list(children = c("CD4_naive", "CD4_Memory"))),
            "NK",
            "Mast"),
        "Endothelial lineage" = c("Endothelial", "Activated endothelial"))

`resolveTreeLeaves()` expands any node, however far up the tree, to the
cell types beneath it. It is used by relative co-localisation when it
conditions on a parent that is broader than a cell type's own group (see
`parentOverrides` in
[`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md)),
and `autoParentOverrides()` finds those parents automatically.

## Examples

``` r
data(lungCellTypeTree)
buildCellTypeMap(lungCellTypeTree)
#> $cellMap
#>                       cellType               group               broad
#> 1         Alveolar Macrophages  Macrophage lineage              Immune
#> 2     Interstitial Macrophages  Macrophage lineage              Immune
#> 3            SPP1+ Macrophages  Macrophage lineage              Immune
#> 4  Macrophages - IFN-activated  Macrophage lineage              Immune
#> 5               Monocytes/MDMs  Macrophage lineage              Immune
#> 6                         cDCs     Dendritic cells              Immune
#> 7                Migratory DCs     Dendritic cells              Immune
#> 8             Langerhans cells     Dendritic cells              Immune
#> 9                         pDCs     Dendritic cells              Immune
#> 10                        Mast        Granulocytes              Immune
#> 11                   Basophils        Granulocytes              Immune
#> 12                 Neutrophils        Granulocytes              Immune
#> 13       Proliferating Myeloid                <NA>              Immune
#> 14                CD4+ T-cells             T cells              Immune
#> 15                CD8+ T-cells             T cells              Immune
#> 16                       Tregs             T cells              Immune
#> 17       Proliferating T-cells             T cells              Immune
#> 18                      NK/NKT          NK lineage              Immune
#> 19        Proliferating NK/NKT          NK lineage              Immune
#> 20                     B cells           B lineage              Immune
#> 21                      Plasma           B lineage              Immune
#> 22       Proliferating B cells           B lineage              Immune
#> 23                Alveolar FBs  Fibroblast lineage         Mesenchymal
#> 24             Adventitial FBs  Fibroblast lineage         Mesenchymal
#> 25      Activated Fibrotic FBs  Fibroblast lineage         Mesenchymal
#> 26            Inflammatory FBs  Fibroblast lineage         Mesenchymal
#> 27              Subpleural FBs  Fibroblast lineage         Mesenchymal
#> 28           Proliferating FBs  Fibroblast lineage         Mesenchymal
#> 29              Myofibroblasts  Fibroblast lineage         Mesenchymal
#> 30              SMCs/Pericytes                <NA>         Mesenchymal
#> 31                 Mesothelial                <NA>         Mesenchymal
#> 32                   Arteriole Endothelial lineage Endothelial lineage
#> 33                      Venous Endothelial lineage Endothelial lineage
#> 34                   Capillary Endothelial lineage Endothelial lineage
#> 35                   Lymphatic Endothelial lineage Endothelial lineage
#> 36                         AT1            Alveolar  Epithelial lineage
#> 37                         AT2            Alveolar  Epithelial lineage
#> 38            Transitional AT2            Alveolar  Epithelial lineage
#> 39           Proliferating AT2            Alveolar  Epithelial lineage
#> 40                       Basal              Airway  Epithelial lineage
#> 41                   Secretory              Airway  Epithelial lineage
#> 42                      Goblet              Airway  Epithelial lineage
#> 43               Multiciliated              Airway  Epithelial lineage
#> 44                        RASC              Airway  Epithelial lineage
#> 45                        PNEC              Airway  Epithelial lineage
#> 46                KRT5-/KRT17+              Airway  Epithelial lineage
#> 47        Proliferating Airway              Airway  Epithelial lineage
#> 
#> $groupsWithChildren
#>  [1] "Macrophage lineage"  "Dendritic cells"     "Granulocytes"       
#>  [4] "T cells"             "NK lineage"          "B lineage"          
#>  [7] "Fibroblast lineage"  "Endothelial lineage" "Alveolar"           
#> [10] "Airway"             
#> 
resolveTreeLeaves(lungCellTypeTree, "Immune")
#>  [1] "Alveolar Macrophages"        "Interstitial Macrophages"   
#>  [3] "SPP1+ Macrophages"           "Macrophages - IFN-activated"
#>  [5] "Monocytes/MDMs"              "cDCs"                       
#>  [7] "Migratory DCs"               "Langerhans cells"           
#>  [9] "pDCs"                        "Mast"                       
#> [11] "Basophils"                   "Neutrophils"                
#> [13] "Proliferating Myeloid"       "CD4+ T-cells"               
#> [15] "CD8+ T-cells"                "Tregs"                      
#> [17] "Proliferating T-cells"       "NK/NKT"                     
#> [19] "Proliferating NK/NKT"        "B cells"                    
#> [21] "Plasma"                      "Proliferating B cells"      
autoParentOverrides(lungCellTypeTree)
#>   Macrophage lineage      Dendritic cells         Granulocytes 
#>            "Myeloid"            "Myeloid"            "Myeloid" 
#>              T cells           NK lineage            B lineage 
#>           "Lymphoid"           "Lymphoid"           "Lymphoid" 
#>   Fibroblast lineage             Alveolar               Airway 
#>        "Mesenchymal" "Epithelial lineage" "Epithelial lineage" 

data(lungSubset)
cells = annotateCellTypeTree(getCellData(lungSubset), lungCellTypeTree)
table(cells$group, useNA = "ifany")
#> 
#>              Airway            Alveolar           B lineage     Dendritic cells 
#>                1899                3438                1210                 193 
#> Endothelial lineage  Fibroblast lineage        Granulocytes  Macrophage lineage 
#>                6214                5148                1575                4571 
#>          NK lineage             T cells                <NA> 
#>                 607                2927                2189 
```
