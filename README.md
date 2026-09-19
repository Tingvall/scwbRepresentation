# scwbRepresentation

This repository is a provisional, explicitly activated scworkbench operation
pack containing bounded Leiden, PCA, UMAP, sparse kNN-neighbor, and kNN
pseudobulk operations.

It is not the complete Representation pack and is not part of a standard
profile. PCA consumes bounded typed matrix slices, authenticates its requested
component axis against matrix rank, and publishes the space plus loadings,
variance, and fit summary as one promotion closure. Installation remains inert;
explicit activation is required. The kNN operations consume separate
authenticated query and candidate space axes. They publish a deterministic
query/rank edge ledger, and pseudobulk reads only bounded feature/candidate
matrix slices tied to that ledger.
