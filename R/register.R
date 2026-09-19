#' Register the provisional scwbRepresentation operation pack
#'
#' Installation is inert. scworkbench invokes this zero-argument entry point
#' only at its explicit operation-pack activation boundary.
#'
#' @return A validated scworkbench operation-pack registration.
#' @export
register_operation_pack <- function() {
  scworkbench::wb_operation_pack_register(
    system.file("scworkbench", "operation-pack.json",
      package = "scwbRepresentation"),
    list(`scwbrepresentation::leiden` = .run_leiden,
      `scwbrepresentation::pca_fit` = .run_pca_fit,
      `scwbrepresentation::umap_fit` = .run_umap_fit,
      `scwbrepresentation::umap_transform` = .run_umap_transform,
      `scwbrepresentation::knn_neighbors` = .run_knn_neighbors,
      `scwbrepresentation::knn_pseudobulk` = .run_knn_pseudobulk))
}
