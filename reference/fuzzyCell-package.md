# fuzzyCell: Bootstrap Classification Confidence for Spatial Cell-Type Analysis

Cell type labels in spatial single-cell data are uncertain, particularly
between closely related subtypes. fuzzyCell estimates a per-cell
classification confidence by bootstrapping (multinomial resampling of
each cell's counts followed by re-classification) and uses the resulting
confidences as soft cell type weights. The weights are propagated
through downstream analyses: cell type proportions, co-localisation
(spicyR's L-function, and Statial's Kontextual test for relative
co-localisation) and niche identification (lisaClust), each compared to
the same analysis using hard labels. Cell types are organised with a
user-supplied hierarchy of lineages, groups and subtypes. Helpers for
testing against sample-level outcomes and for plotting are included.

## See also

Useful links:

- <https://github.com/shreyarajeshrao/fuzzyCell>

- Report bugs at <https://github.com/shreyarajeshrao/fuzzyCell/issues>

## Author

**Maintainer**: Shreya Rao <shreya.rajeshrao@sydney.edu.au>
([ORCID](https://orcid.org/0009-0007-8140-6844))

Authors:

- Shreya Rao <shreya.rajeshrao@sydney.edu.au>
  ([ORCID](https://orcid.org/0009-0007-8140-6844))

- Ellis Patrick \[thesis advisor\]
