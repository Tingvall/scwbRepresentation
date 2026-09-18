representation_activate_pack <- function() {
  library(scworkbench)
  if ("scwbrepresentation::leiden" %in% wb_list_ops()) return(invisible(TRUE))
  path <- system.file("scworkbench", "operation-pack.json",
    package = "scwbRepresentation")
  if (!nzchar(path)) stop("scwbRepresentation is not installed")
  wb_operation_pack_activate(path, lib.loc = .libPaths())
  invisible(TRUE)
}

representation_fixture <- function() {
  representation_activate_pack()
  root <- tempfile("scwbRepresentation-test-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  cells <- paste0("c", seq_len(8L))
  coordinates <- cbind(x = c(0, .1, .2, .1, 10, 10.1, 10.2, 10.1),
    y = c(0, .1, 0, -.1, 10, 10.1, 10, 9.9),
    z = c(1, 1.1, .9, 1, 3, 3.1, 2.9, 3))
  rownames(coordinates) <- cells
  wb_define_cells(bundle, cells)
  wb_add_space(bundle, "latent", coordinates, overwrite = TRUE)
  space <- structure(list(kind = "space", space_id = "latent",
    values = coordinates, dimensions = as.integer(dim(coordinates)),
    cell_id = cells, storage_profile = "space_dense_v1"),
    class = "wb_primitive")
  list(root = root, bundle = bundle, cells = cells, space = space)
}

representation_execute <- function(fixture, inputs = list()) {
  transport <- getFromNamespace(".wb_operation_pack_params_transport",
    "scworkbench")
  getFromNamespace(".wb_operation_pack_v2_execute", "scworkbench")(
    fixture$bundle, "scwbrepresentation::leiden",
    list(space = fixture$space), transport(inputs),
    effective_cells = fixture$cells)
}
