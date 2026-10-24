# scwbRepresentation

This repository is an explicitly activated scworkbench ABI-v3 operation
pack containing bounded Leiden, PCA, UMAP, sparse kNN-neighbor, and kNN
pseudobulk operations.

The standard profile selects this bounded six-operation pack. It is not a claim
to implement every representation method. PCA consumes bounded typed matrix slices, authenticates its requested
component axis against matrix rank, and publishes the space plus loadings,
variance, and fit summary as one promotion closure. Installation remains inert;
explicit activation is required. The kNN operations consume separate
authenticated query and candidate space axes. They publish a deterministic
query/rank edge ledger, and pseudobulk reads only bounded feature/candidate
matrix slices tied to that ledger.

```r
library(scworkbench)
wb_operation_pack_load("scwbRepresentation")
wb_op_interface("scwbrepresentation::pca_fit")
```

The six operations are `leiden`, `pca_fit`, `umap_fit`, `umap_transform`,
`knn_neighbors`, and `knn_pseudobulk`, all under `scwbrepresentation::`.
Inspect the interface for parameters, typed references and terminal closure.
Core supplies slicing and owns publication/provenance; pack code receives no
committed storage authority. Inspection is pack-free; execution/replay requires
the compatible activated pack.

## PCA variance denominators

`variance` and `fit_summary` retain the historical fields unchanged.
`variance_explained` divides squared component standard deviations by the sum
of **all sdevs returned by the fit backend**, identified by
`variance_explained_denominator` and its numeric `_value`. For truncated irlba
that is the retained block; prcomp can return a full spectrum even when scores
are truncated. Do not interpret every historical cumulative value of one as
all input variation captured.

`retained_component_variance_fraction` and its `cumulative_` counterpart divide
by `retained_component_variance_total` (only emitted components). Separately,
`input_total_variance_fraction` and `cumulative_input_total_variance_fraction`
divide by `input_total_variance`: the sum of squared entries of the matrix
actually fitted, after constant-feature removal, optional centering and backend
scaling, divided by `n_cells - 1`. This is total sample variance when centered;
without centering it is a second-moment quantity including the mean. With
scaling it equals the number of fitted nonconstant features for both supported
backends. Existing feature moments supply this total without another input scan.
A nonfinite/nonpositive denominator yields numeric missing input-total fractions
and `input_total_variance_status="nonfinite_or_nonpositive_input_total"`.

For a two-component truncated fit, retained cumulative fraction is one while
input-total cumulative fraction can be much smaller. Inspect the variance
obs-table using the promoted `variance` artifact and `wb_art_load()$table`.
New outputs have new content/pack identities; old artifacts and their identity
schemes are not rewritten. This changed manifest requires its own qualification;
the old standard-profile receipt does not qualify this source tuple.
