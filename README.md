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
