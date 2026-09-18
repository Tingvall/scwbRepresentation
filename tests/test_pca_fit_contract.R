test_representation_pca_fit_parameter_axis_contract <- function() {
  representation_activate_pack()
  root <- tempfile("scwbRepresentation-pca-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  on.exit({try(wb_close_bundle(bundle), silent = TRUE);
    unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
  cells <- paste0("c", seq_len(12L)); features <- paste0("g", seq_len(5L))
  set.seed(17L)
  values <- matrix(sample.int(40L, length(cells) * length(features),
    replace = TRUE), nrow = length(features),
    dimnames = list(features, cells))
  wb_define_cells(bundle, cells)
  wb_attach_matrix(bundle, "RNA", values, feature_type = "gene")
  run <- wb_session_run_op(wb_session_new(), bundle, "pca_fit",
    inputs = list(matrix = "RNA", cell_set = "all", feature_ids = features,
      n_components = 3L, center = TRUE, scale = FALSE,
      backend = "prcomp", seed = 1L, space_id_out = "pca",
      overwrite = FALSE, max_dense_gb = 1, chunk_size = 2L))
  output <- run$run$outputs
  stopifnot(identical(dim(output$space$values), c(12L, 3L)),
    identical(rownames(output$space$values), cells),
    identical(colnames(output$space$values), paste0("PC", 1:3)),
    nrow(output$loadings$table) == 15L,
    identical(output$loadings$table$feature_id, rep(features, 3L)),
    identical(output$loadings$table$PC, rep(paste0("PC", 1:3),
      each = 5L)), nrow(output$variance$table) == 3L,
    output$fit_summary$table$n_components[[1]] == 3L)
  promoted <- wb_session_promote_current(bundle, run$session,
    outputs = "fit_summary")
  stopifnot(all(promoted$promoted$status == "ok"),
    identical(dim(wb_space(bundle, "pca")), c(12L, 3L)))
  before_replace <- wb_space(bundle, "pca")

  execute_bad <- function(mutator) {
    options(scworkbench.operation_pack_v3_output_test_hook = mutator)
    on.exit(options(scworkbench.operation_pack_v3_output_test_hook = NULL),
      add = TRUE)
    try(wb_session_run_op(wb_session_new(), bundle, "pca_fit",
      inputs = list(matrix = "RNA", cell_set = "all",
        feature_ids = features, n_components = 3L, backend = "prcomp",
        scale = FALSE, space_id_out = "bad", max_dense_gb = 1)),
      silent = TRUE)
  }
  truncate_space <- function(n) function(writes, ...) {
    writes$space$value$values <- writes$space$value$values[, seq_len(n),
      drop = FALSE]
    writes
  }
  four_space <- function(writes, ...) {
    values <- writes$space$value$values
    writes$space$value$values <- cbind(values, PC4 = values[, 1L])
    writes
  }
  missing_cell <- function(writes, ...) {
    writes$space$value$values <- writes$space$value$values[-1L, , drop = FALSE]
    writes
  }
  duplicate_dimension <- function(writes, ...) {
    colnames(writes$space$value$values)[[2L]] <- "PC1"
    writes
  }
  partial <- function(writes, ...) { writes$variance <- NULL; writes }
  stopifnot(inherits(execute_bad(truncate_space(2L)), "try-error"),
    inherits(execute_bad(four_space), "try-error"),
    inherits(execute_bad(missing_cell), "try-error"),
    inherits(execute_bad(duplicate_dimension), "try-error"),
    inherits(execute_bad(partial), "try-error"),
    !"bad" %in% wb_spaces(bundle)$space_id)

  options(scworkbench.operation_pack_v3_output_test_hook = partial)
  failed_replace <- try(wb_session_run_op(wb_session_new(), bundle,
    "pca_fit", inputs = list(matrix = "RNA", cell_set = "all",
      feature_ids = features, n_components = 3L, center = FALSE,
      scale = FALSE, backend = "prcomp", space_id_out = "pca",
      overwrite = TRUE, max_dense_gb = 1)), silent = TRUE)
  options(scworkbench.operation_pack_v3_output_test_hook = NULL)
  stopifnot(inherits(failed_replace, "try-error"),
    identical(before_replace, wb_space(bundle, "pca")))
  replacement <- wb_session_run_op(wb_session_new(), bundle, "pca_fit",
    inputs = list(matrix = "RNA", cell_set = "all",
      feature_ids = features, n_components = 3L, center = FALSE,
      scale = FALSE, backend = "prcomp", space_id_out = "pca",
      overwrite = TRUE, max_dense_gb = 1))
  replacement_promotion <- wb_session_promote_current(bundle,
    replacement$session, outputs = "fit_summary")
  stopifnot(all(replacement_promotion$promoted$status == "ok"),
    !identical(before_replace, wb_space(bundle, "pca")),
    isTRUE(wb_bundle_validation_report(root, "content")$ok))
  TRUE
}

test_representation_pca_fit_old_core_equivalence_contract <- function() {
  representation_activate_pack()
  fixture_path <- if (file.exists(file.path("tests", "fixtures",
      "pca_old_core_reference.json"))) file.path("tests", "fixtures",
      "pca_old_core_reference.json") else file.path("fixtures",
      "pca_old_core_reference.json")
  reference <- jsonlite::read_json(fixture_path, simplifyVector = FALSE)
  root <- tempfile("scwbRepresentation-pca-equivalence-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  on.exit({try(wb_close_bundle(bundle), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
  cells <- as.character(unlist(reference$cells, use.names = FALSE))
  features <- as.character(unlist(reference$features, use.names = FALSE))
  values <- do.call(rbind, lapply(reference$values, function(x)
    as.numeric(unlist(x, use.names = FALSE))))
  dimnames(values) <- list(features, cells)
  wb_define_cells(bundle, cells)
  wb_attach_matrix(bundle, "RNA", values, feature_type = "gene")
  result <- wb_session_run_op(wb_session_new(), bundle, "pca_fit",
    inputs = list(matrix = "RNA", cell_set = "all", feature_ids = features,
      n_components = 3L, center = TRUE, scale = FALSE,
      backend = "prcomp", seed = 1L, space_id_out = "pca",
      overwrite = FALSE, max_dense_gb = 1, chunk_size = 100L))
  view <- wb_session_view(result$session)
  observed <- lapply(result$run$outputs, function(x) {
    identity <- scworkbench:::.wb_prov_payload_identity(view, x)
    list(scheme = identity$content_identity_scheme,
      sha256 = identity$content_sha256)
  })
  stopifnot(identical(observed[names(reference$outputs)],
      reference$outputs),
    identical(result$run$params$n_components, 3L),
    identical(as.character(unlist(result$run$params$feature_ids,
      use.names = FALSE)), features),
    identical(result$run$params$space_id_out, "pca"),
    identical(result$run$params$overwrite, FALSE))
  TRUE
}
