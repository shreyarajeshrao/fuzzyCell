# Niche plots

`plotNicheSpatial()` shows the cells of one niche in one image.
`nicheCompositionMap()` shows the relative enrichment of each cell type
in each niche as a bubble plot or heatmap. `nicheBarplot()` shows it for
a single niche. `nicheProportionBoxplot()` compares a niche's per-image
proportion between outcome groups.

## Usage

``` r
plotNicheSpatial(spe, imageID, niche, nicheCol = "nicheWeighted", cellSize = 1)

nicheCompositionMap(
  spe,
  niche,
  cellType = "cellType",
  type = c("bubble", "heatmap"),
  limit = c(0.33, 3),
  ...
)

nicheBarplot(spe, niche, target, cellType = "cellType", limit = c(0.33, 3))

nicheProportionBoxplot(proportions, niche, spe, outcomeCol)
```

## Arguments

- spe:

  A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with the cell
  type, image and coordinates of every cell in the form given by
  [`getCellData()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/getCellData.md).

- imageID:

  The image to plot.

- niche:

  For `plotNicheSpatial()` and the proportion boxplot, the niche to
  show. For `nicheCompositionMap()` and `nicheBarplot()`, the niche
  assignments: either a vector of niches named by cell ID, such as
  `niches$nicheWeighted` from
  [`runNicheIdentification()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nicheIdentification.md),
  or the name of a `colData` column holding them.

- nicheCol:

  The niche assignments for `plotNicheSpatial()`, in either of the forms
  taken by `niche` in `nicheCompositionMap()`.

- cellSize:

  Size of the highlighted cells in `plotNicheSpatial()`, the cells of
  the chosen niche. The grey background cells always have size 0.7.

- cellType:

  Column of `colData(spe)` with cell type labels.

- type:

  `"bubble"` or `"heatmap"`.

- limit:

  Lower and upper limit of the relative enrichment.

- ...:

  Further arguments passed to
  [`pheatmap::pheatmap()`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html).

- target:

  The niche to plot in `nicheBarplot()`.

- proportions:

  The `propWeighted` or `propUnweighted` element of the results of
  [`testNicheOutcome()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/testNicheOutcome.md):
  the table you pass decides whether the weighted or the unweighted
  niches are plotted.

- outcomeCol:

  Name of the per-image outcome column in `colData(spe)`.

## Value

A `ggplot` object; for `nicheProportionBoxplot()` a `girafe`
(interactive) widget; for `nicheCompositionMap(type = "heatmap")` the
result of
[`pheatmap::pheatmap()`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html).

## Examples

``` r
data(lungSubset)
data(lungCellTypeTree)
data(lungConfList)
niches = runNicheIdentification(lungSubset, lungCellTypeTree, lungConfList,
    r = 50, k = 6)
#> Computing unweighted LISA curves...
#> Generating local indicators of spatial association.
#> Computing weighted LISA curves...
#> Generating local indicators of spatial association.
#> Clustering niches...
plotNicheSpatial(lungSubset, "VUHD113", "niche1",
    nicheCol = niches$nicheWeighted)

nicheCompositionMap(lungSubset, niches$nicheWeighted)

nicheBarplot(lungSubset, niches$nicheWeighted, "niche1")

res = testNicheOutcome(lungSubset, niches$nicheUnweighted,
    niches$nicheWeighted, outcomeCol = "disease_status", subject = "subject")
if (requireNamespace("ggiraph", quietly = TRUE)) {
    nicheProportionBoxplot(res$propWeighted, "niche1", lungSubset,
        "disease_status")
    nicheProportionBoxplot(res$propUnweighted, "niche1", lungSubset,
        "disease_status")
}

{"x":{"html":"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' class='ggiraph-svg' role='graphics-document' id='svg_e5c8c404fe174e4c' viewBox='0 0 432 252'>\n <defs id='svg_e5c8c404fe174e4c_defs'>\n  <clipPath id='svg_e5c8c404fe174e4c_c1'>\n   <rect x='0' y='0' width='432' height='252'/>\n  <\/clipPath>\n  <clipPath id='svg_e5c8c404fe174e4c_c2'>\n   <rect x='35.24' y='22.77' width='314.82' height='210.94'/>\n  <\/clipPath>\n <\/defs>\n <g id='svg_e5c8c404fe174e4c_rootg' class='ggiraph-svg-rootg'>\n  <g clip-path='url(#svg_e5c8c404fe174e4c_c1)'>\n   <rect x='0' y='0' width='432' height='252' fill='#FFFFFF' fill-opacity='1' stroke='#FFFFFF' stroke-opacity='1' stroke-width='0.75' stroke-linejoin='round' stroke-linecap='round' class='ggiraph-svg-bg'/>\n   <rect x='0' y='0' width='432' height='252' fill='#FFFFFF' fill-opacity='1' stroke='#FFFFFF' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='round'/>\n  <\/g>\n  <g clip-path='url(#svg_e5c8c404fe174e4c_c2)'>\n   <rect x='35.24' y='22.77' width='314.82' height='210.94' fill='#FFFFFF' fill-opacity='1' stroke='none'/>\n   <line x1='121.1' y1='118.26' x2='121.1' y2='88.11' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <line x1='121.1' y1='195.1' x2='121.1' y2='224.11' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <polygon points='85.32,118.26 85.32,195.10 156.87,195.10 156.87,118.26 85.32,118.26' fill='#FFFFFF' fill-opacity='1' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <line x1='85.32' y1='190.55' x2='156.87' y2='190.55' stroke='#F8766D' stroke-opacity='1' stroke-width='2.13' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <line x1='264.2' y1='108.17' x2='264.2' y2='108.17' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <line x1='264.2' y1='113.22' x2='264.2' y2='113.22' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <polygon points='228.42,108.17 228.42,113.22 299.97,113.22 299.97,108.17 228.42,108.17' fill='#FFFFFF' fill-opacity='1' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <line x1='228.42' y1='109.82' x2='299.97' y2='109.82' stroke='#00BFC4' stroke-opacity='1' stroke-width='2.13' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <circle id='svg_e5c8c404fe174e4c_e1' cx='249.99' cy='206.01' r='1.47pt' fill='#00BFC4' fill-opacity='0.7' stroke='#00BFC4' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='TILD049MA' data-id='TILD049MA'/>\n   <circle id='svg_e5c8c404fe174e4c_e2' cx='250.03' cy='109.4' r='1.47pt' fill='#00BFC4' fill-opacity='0.7' stroke='#00BFC4' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='TILD113LA' data-id='TILD113LA'/>\n   <circle id='svg_e5c8c404fe174e4c_e3' cx='271.75' cy='108.66' r='1.47pt' fill='#00BFC4' fill-opacity='0.7' stroke='#00BFC4' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='TILD167LA' data-id='TILD167LA'/>\n   <circle id='svg_e5c8c404fe174e4c_e4' cx='112.59' cy='189.93' r='1.47pt' fill='#F8766D' fill-opacity='0.7' stroke='#F8766D' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='VUHD038' data-id='VUHD038'/>\n   <circle id='svg_e5c8c404fe174e4c_e5' cx='107' cy='118.56' r='1.47pt' fill='#F8766D' fill-opacity='0.7' stroke='#F8766D' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='VUHD069' data-id='VUHD069'/>\n   <circle id='svg_e5c8c404fe174e4c_e6' cx='107.51' cy='224.12' r='1.47pt' fill='#F8766D' fill-opacity='0.7' stroke='#F8766D' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='VUHD090' data-id='VUHD090'/>\n   <circle id='svg_e5c8c404fe174e4c_e7' cx='116.14' cy='88.71' r='1.47pt' fill='#F8766D' fill-opacity='0.7' stroke='#F8766D' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='VUHD095' data-id='VUHD095'/>\n   <circle id='svg_e5c8c404fe174e4c_e8' cx='130.68' cy='194.61' r='1.47pt' fill='#F8766D' fill-opacity='0.7' stroke='#F8766D' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='VUHD113' data-id='VUHD113'/>\n   <circle id='svg_e5c8c404fe174e4c_e9' cx='267.87' cy='113.02' r='1.47pt' fill='#00BFC4' fill-opacity='0.7' stroke='#00BFC4' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='VUILD105MA2' data-id='VUILD105MA2'/>\n   <circle id='svg_e5c8c404fe174e4c_e10' cx='277.43' cy='32.36' r='1.47pt' fill='#00BFC4' fill-opacity='0.7' stroke='#00BFC4' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round' title='VUILD142MA' data-id='VUILD142MA'/>\n  <\/g>\n  <g clip-path='url(#svg_e5c8c404fe174e4c_c1)'>\n   <polyline points='35.24,233.71 35.24,22.77' fill='none' stroke='#000000' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='square'/>\n   <text x='18.07' y='175.9' font-size='6.6pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>0.4<\/text>\n   <text x='18.07' y='114.18' font-size='6.6pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>0.5<\/text>\n   <text x='18.07' y='52.46' font-size='6.6pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>0.6<\/text>\n   <polyline points='32.50,172.88 35.24,172.88' fill='none' stroke='#000000' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <polyline points='32.50,111.16 35.24,111.16' fill='none' stroke='#000000' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <polyline points='32.50,49.44 35.24,49.44' fill='none' stroke='#000000' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <polyline points='35.24,233.71 350.06,233.71' fill='none' stroke='#000000' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='square'/>\n   <polyline points='121.10,236.44 121.10,233.71' fill='none' stroke='#000000' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <polyline points='264.20,236.44 264.20,233.71' fill='none' stroke='#000000' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='round' stroke-linecap='butt'/>\n   <text x='106.91' y='244.69' font-size='6.6pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>Control<\/text>\n   <text x='248.3' y='244.69' font-size='6.6pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>Disease<\/text>\n   <text transform='translate(13.05,170.75) rotate(-90.00)' font-size='8.25pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>niche1 proportion<\/text>\n   <rect x='361.01' y='105.48' width='65.51' height='45.52' fill='#FFFFFF' fill-opacity='1' stroke='none'/>\n   <rect x='366.49' y='110.96' width='17.28' height='17.28' fill='#FFFFFF' fill-opacity='1' stroke='none'/>\n   <polyline points='375.13,126.51 375.13,123.92' fill='none' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='375.13,115.28 375.13,112.68' fill='none' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <rect x='368.65' y='115.28' width='12.96' height='8.64' fill='#FFFFFF' fill-opacity='1' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='368.65,119.60 381.61,119.60' fill='none' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='375.13,126.51 375.13,126.51' fill='none' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='375.13,112.68 375.13,112.68' fill='none' stroke='#F8766D' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <circle cx='375.13' cy='119.6' r='1.47pt' fill='#F8766D' fill-opacity='0.7' stroke='#F8766D' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round'/>\n   <rect x='366.49' y='128.24' width='17.28' height='17.28' fill='#FFFFFF' fill-opacity='1' stroke='none'/>\n   <polyline points='375.13,143.79 375.13,141.20' fill='none' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='375.13,132.56 375.13,129.96' fill='none' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <rect x='368.65' y='132.56' width='12.96' height='8.64' fill='#FFFFFF' fill-opacity='1' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='368.65,136.88 381.61,136.88' fill='none' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='375.13,143.79 375.13,143.79' fill='none' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <polyline points='375.13,129.96 375.13,129.96' fill='none' stroke='#00BFC4' stroke-opacity='1' stroke-width='1.07' stroke-linejoin='miter' stroke-linecap='butt'/>\n   <circle cx='375.13' cy='136.88' r='1.47pt' fill='#00BFC4' fill-opacity='0.7' stroke='#00BFC4' stroke-opacity='0.7' stroke-width='0.71' stroke-linejoin='round' stroke-linecap='round'/>\n   <text x='389.25' y='122.62' font-size='6.6pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>Control<\/text>\n   <text x='389.25' y='139.9' font-size='6.6pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>Disease<\/text>\n   <text x='35.24' y='14.55' font-size='9.9pt' font-family='Liberation Sans, Arial, Helvetica, sans-serif'>niche1 (unweighted)<\/text>\n  <\/g>\n <\/g>\n<\/svg>","js":null,"uid":"svg_e5c8c404fe174e4c","ratio":1.714285714285714,"settings":{"tooltip":{"css":".tooltip_SVGID_ { padding:5px;background:black;color:white;border-radius:2px;text-align:left; ; position:absolute;pointer-events:none;z-index:9999;}","placement":"doc","opacity":0.9,"offx":10,"offy":10,"use_cursor_pos":true,"use_fill":false,"use_stroke":false,"delay_over":200,"delay_out":500},"hover":{"css":".hover_data_SVGID_ { fill:orange;stroke:black;cursor:pointer; }\ntext.hover_data_SVGID_ { stroke:none;fill:orange; }\ncircle.hover_data_SVGID_ { fill:orange;stroke:black; }\nline.hover_data_SVGID_, polyline.hover_data_SVGID_ { fill:none;stroke:orange; }\nrect.hover_data_SVGID_, polygon.hover_data_SVGID_, path.hover_data_SVGID_ { fill:orange;stroke:none; }\nimage.hover_data_SVGID_ { stroke:orange; }","reactive":true,"nearest_distance":null,"linked":false},"hover_inv":{"css":""},"hover_key":{"css":".hover_key_SVGID_ { fill:orange;stroke:black;cursor:pointer; }\ntext.hover_key_SVGID_ { stroke:none;fill:orange; }\ncircle.hover_key_SVGID_ { fill:orange;stroke:black; }\nline.hover_key_SVGID_, polyline.hover_key_SVGID_ { fill:none;stroke:orange; }\nrect.hover_key_SVGID_, polygon.hover_key_SVGID_, path.hover_key_SVGID_ { fill:orange;stroke:none; }\nimage.hover_key_SVGID_ { stroke:orange; }","reactive":true},"hover_theme":{"css":".hover_theme_SVGID_ { fill:orange;stroke:black;cursor:pointer; }\ntext.hover_theme_SVGID_ { stroke:none;fill:orange; }\ncircle.hover_theme_SVGID_ { fill:orange;stroke:black; }\nline.hover_theme_SVGID_, polyline.hover_theme_SVGID_ { fill:none;stroke:orange; }\nrect.hover_theme_SVGID_, polygon.hover_theme_SVGID_, path.hover_theme_SVGID_ { fill:orange;stroke:none; }\nimage.hover_theme_SVGID_ { stroke:orange; }","reactive":true},"select":{"css":".select_data_SVGID_ { fill:red;stroke:black;cursor:pointer; }\ntext.select_data_SVGID_ { stroke:none;fill:red; }\ncircle.select_data_SVGID_ { fill:red;stroke:black; }\nline.select_data_SVGID_, polyline.select_data_SVGID_ { fill:none;stroke:red; }\nrect.select_data_SVGID_, polygon.select_data_SVGID_, path.select_data_SVGID_ { fill:red;stroke:none; }\nimage.select_data_SVGID_ { stroke:red; }","type":"multiple","only_shiny":true,"selected":[],"linked":false},"select_inv":{"css":""},"select_key":{"css":".select_key_SVGID_ { fill:red;stroke:black;cursor:pointer; }\ntext.select_key_SVGID_ { stroke:none;fill:red; }\ncircle.select_key_SVGID_ { fill:red;stroke:black; }\nline.select_key_SVGID_, polyline.select_key_SVGID_ { fill:none;stroke:red; }\nrect.select_key_SVGID_, polygon.select_key_SVGID_, path.select_key_SVGID_ { fill:red;stroke:none; }\nimage.select_key_SVGID_ { stroke:red; }","type":"single","only_shiny":true,"selected":[]},"select_theme":{"css":".select_theme_SVGID_ { fill:red;stroke:black;cursor:pointer; }\ntext.select_theme_SVGID_ { stroke:none;fill:red; }\ncircle.select_theme_SVGID_ { fill:red;stroke:black; }\nline.select_theme_SVGID_, polyline.select_theme_SVGID_ { fill:none;stroke:red; }\nrect.select_theme_SVGID_, polygon.select_theme_SVGID_, path.select_theme_SVGID_ { fill:red;stroke:none; }\nimage.select_theme_SVGID_ { stroke:red; }","type":"single","only_shiny":true,"selected":[]},"zoom":{"min":1,"max":1,"duration":300,"default_on":false},"toolbar":{"position":"topright","pngname":"diagram","tooltips":null,"fixed":false,"hidden":[],"delay_over":200,"delay_out":500},"sizing":{"rescale":true,"width":1}}},"evals":[],"jsHooks":[]}
```
