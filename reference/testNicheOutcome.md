# Test niche proportions against an outcome

Computes the proportion of each image's cells in each niche and tests
every niche against a two-level outcome with a Wilcoxon or t-test.

## Usage

``` r
testNicheOutcome(
  spe,
  nicheUnweighted,
  nicheWeighted,
  outcomeCol,
  test = c("wilcox", "ttest"),
  subject = NULL,
  cellID = "cellID",
  imageID = "imageID"
)
```

## Arguments

- spe:

  A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with a cell ID,
  an image and the outcome for every cell in `colData`.

- nicheUnweighted, nicheWeighted:

  Named cell-to-niche vectors, as returned by
  [`runNicheIdentification()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nicheIdentification.md).

- outcomeCol:

  Name of the per-image outcome column in `colData(spe)`. It must have
  exactly two levels (subset your data to two conditions first) and a
  single value per image: repeated measures and mixed models are not
  supported. An error is raised otherwise, before any computation.

- test:

  `"wilcox"` or `"ttest"`.

- subject:

  Name of the `colData(spe)` column identifying the subject (for example
  the patient) each image comes from, or `NULL`. If given, every image
  must belong to one subject and every subject must contribute only one
  image, as in
  [`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md).
  With `NULL` this cannot be checked, and a message says so.

- cellID, imageID:

  Names of the `colData(spe)` columns holding a unique cell ID (matching
  the names of the niche vectors) and the image each cell is in.

## Value

A list with `testUnweighted`, `testWeighted` (one row per niche) and
`propUnweighted`, `propWeighted` (per image niche proportions, to pass
to
[`nicheProportionBoxplot()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/nichePlots.md)).

## Details

The `estimate` is the difference in the niche proportion for the second
level of the outcome compared with the first (alphabetical order for
text, factor order for factors), for both tests, as in
[`runColocalisationModels()`](https://shreyarajeshrao.github.io/fuzzyCell/reference/runColocalisationModels.md).
For `test = "wilcox"` this is the Hodges-Lehmann estimate of the shift
in location. Each image is one observation, so every image must come
from a different subject if `subject` is given.

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
res = testNicheOutcome(lungSubset, niches$nicheUnweighted,
    niches$nicheWeighted, outcomeCol = "disease_status", subject = "subject")
res$testWeighted
#>    niche    estimate statistic     p.value      p.adj
#> 1 niche4 -0.35213399        25 0.007936508 0.04761905
#> 2 niche6  0.13213213         5 0.150793651 0.37142857
#> 3 niche2  0.06735556         6 0.222222222 0.37142857
#> 4 niche1  0.10960504         7 0.309523810 0.37142857
#> 5 niche5  0.01701702         7 0.309523810 0.37142857
#> 6 niche3  0.00000000        12 1.000000000 1.00000000
```
