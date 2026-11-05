# scwbRepresentation 0.1.1

- Makes bounded-query exact RANN search the default backend for
  `scwbrepresentation::knn_neighbors`.
- Retains the brute-force backend as an explicit fallback and numerical oracle.
- Canonicalizes neighbor ties by authenticated candidate order and falls back
  to the oracle for ambiguous boundary ties.
- Records backend and query chunk size in run provenance without changing the
  neighbor terminal schemas.
- Leaves `scwbrepresentation::knn_pseudobulk` on its established brute-force
  scientific implementation.
