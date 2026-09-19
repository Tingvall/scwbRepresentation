#!/usr/bin/env Rscript

library(scworkbench)
root <- normalizePath(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[
  grepl("^--file=", commandArgs(FALSE))])), ".."), winslash = "/",
  mustWork = TRUE)
bundle_root <- tempfile("scwb-umap-old-core-")
bundle <- wb_init_bundle(bundle_root, overwrite = TRUE)
on.exit({ try(wb_close_bundle(bundle), silent = TRUE)
  unlink(bundle_root, recursive = TRUE, force = TRUE) }, add = TRUE)
cells <- paste0("c", sprintf("%02d", seq_len(30L)))
values <- outer(seq_along(cells), seq_len(4L), function(i, j)
  sin(i * j / 7) + cos(i / (j + 1)))
dimnames(values) <- list(cells, paste0("PC", seq_len(4L)))
wb_define_cells(bundle, cells)
wb_add_space(bundle, "pca", values)
result <- wb_session_run_op(wb_session_new(), bundle, "umap_fit",
  inputs = list(space = "pca", dims = 1:3, space_id_out = "umap.pca",
    n_neighbors = 5L, min_dist = .2, metric = "euclidean",
    nn_method = "nndescent", spread = 1, seed = 17L))
coordinates <- wb_space(wb_session_view(result$session), "umap.pca")
fixture <- list(schema = "scwbRepresentation.umap_old_core_reference.v1",
  old_core_commit = "651cca3cb0f4b1cd18281cbbee982f267b265e6c",
  cells = unname(as.list(cells)), dimensions = unname(as.list(colnames(values))),
  source = lapply(seq_len(nrow(values)), function(i)
    unname(as.list(values[i, ]))),
  coordinates = lapply(seq_len(nrow(coordinates)), function(i)
    unname(as.list(coordinates[i, c("x", "y")]))),
  retained_effective_parameters = result$run$params)
jsonlite::write_json(fixture, file.path(root, "tests", "fixtures",
  "umap_old_core_reference.json"), auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = NA)
