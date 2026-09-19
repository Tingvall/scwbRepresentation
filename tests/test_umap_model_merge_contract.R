representation_umap_fixture <- function() {
  representation_activate_pack()
  root <- tempfile("scwbRepresentation-umap-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  cells <- paste0("c", sprintf("%02d", seq_len(30L)))
  values <- outer(seq_along(cells), seq_len(4L), function(i, j)
    sin(i * j / 7) + cos(i / (j + 1)))
  dimnames(values) <- list(cells, paste0("PC", seq_len(4L)))
  wb_define_cells(bundle, cells)
  wb_add_space(bundle, "pca", values)
  wb_cell_set_replace(bundle, "set:query", cells[seq_len(8L)],
    source = "phase4d_test")
  list(root = root, bundle = bundle, cells = cells, values = values)
}

test_representation_umap_manifest_binding_negatives_contract <- function() {
  representation_activate_pack()
  path <- system.file("scworkbench", "operation-pack.json",
    package = "scwbRepresentation")
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  mutate_and_validate <- function(change) {
    candidate <- change(unserialize(serialize(manifest, NULL)))
    output <- tempfile(fileext = ".json")
    on.exit(unlink(output, force = TRUE), add = TRUE)
    jsonlite::write_json(candidate, output, auto_unbox = TRUE,
      null = "null", digits = NA)
    tryCatch(wb_operation_pack_manifest_validate(output), error = identity)
  }
  fit_input <- mutate_and_validate(function(x) {
    x$operations[[3L]]$publication_policies[[1L]]$model_binding <-
      list(source = "input", name = "space")
    x
  })
  transform_sibling <- mutate_and_validate(function(x) {
    x$operations[[4L]]$publication_policies[[1L]]$model_binding <-
      list(source = "sibling_output", name = "space_2d")
    x
  })
  transform_add <- mutate_and_validate(function(x) {
    x$operations[[4L]]$publication_policies[[1L]]$mode <-
      list(policy = "add_only_v1")
    x
  })
  stopifnot(all(vapply(list(fit_input, transform_sibling, transform_add),
    function(x) inherits(x, "wb_operation_pack_v3_error") && identical(
      x$code, "OPERATION_PACK_MANIFEST_SEMANTICS_INVALID"), logical(1))))
  TRUE
}

test_representation_umap_output_adversarial_contract <- function() {
  if (!requireNamespace("uwot", quietly = TRUE) ||
      !requireNamespace("rnndescent", quietly = TRUE))
    return(structure(TRUE,
      skip = "uwot and rnndescent are required for UMAP adversarial tests"))
  fixture <- representation_umap_fixture(); b <- fixture$bundle
  on.exit({ options(scworkbench.operation_pack_v3_output_test_hook = NULL)
    try(wb_close_bundle(b), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE) }, add = TRUE)
  run_bad <- function(target, hook) {
    options(scworkbench.operation_pack_v3_output_test_hook = hook)
    on.exit(options(scworkbench.operation_pack_v3_output_test_hook = NULL))
    error <- try(wb_session_run_op(wb_session_new(), b, "umap_fit",
      inputs = list(space = "pca", dims = 1:3, space_id_out = target,
        n_neighbors = 5L, nn_method = "nndescent", seed = 23L)),
      silent = TRUE)
    stopifnot(inherits(error, "try-error"),
      !target %in% wb_spaces(b)$space_id)
  }
  run_bad("bad.missing.cell", function(writes, ...) {
    writes$space_2d$value$values <- writes$space_2d$value$values[-1L, ,
      drop = FALSE]; writes
  })
  run_bad("bad.foreign.cell", function(writes, ...) {
    rownames(writes$space_2d$value$values)[[1L]] <- "foreign"
    writes$space_2d$value$cell_id[[1L]] <- "foreign"; writes
  })
  run_bad("bad.duplicate.cell", function(writes, ...) {
    rownames(writes$space_2d$value$values)[[2L]] <-
      rownames(writes$space_2d$value$values)[[1L]]
    writes$space_2d$value$cell_id[[2L]] <-
      writes$space_2d$value$cell_id[[1L]]; writes
  })
  run_bad("bad.dimensions", function(writes, ...) {
    writes$space_2d$value$values <- cbind(writes$space_2d$value$values,
      z = 1); writes
  })
  run_bad("bad.model.extra", function(writes, ...) {
    writeBin(as.raw(1:4), file.path(writes$model$tree_root, "extra.bin"))
    writes
  })
  run_bad("bad.model.missing", function(writes, ...) {
    unlink(file.path(writes$model$tree_root,
      writes$model$tree_manifest$primary_path), force = TRUE)
    writes
  })
  run_bad("bad.model.changed", function(writes, ...) {
    path <- file.path(writes$model$tree_root,
      writes$model$tree_manifest$primary_path)
    writeBin(as.raw(0L), path, useBytes = TRUE)
    writes
  })
  TRUE
}

test_representation_umap_stale_preimage_and_retry_contract <- function() {
  if (!requireNamespace("uwot", quietly = TRUE) ||
      !requireNamespace("rnndescent", quietly = TRUE))
    return(structure(TRUE,
      skip = "uwot and rnndescent are required for the stale-preimage test"))
  fixture <- representation_umap_fixture(); b <- fixture$bundle
  on.exit({ try(wb_close_bundle(b), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE) }, add = TRUE)
  fit <- wb_session_run_op(wb_session_new(), b, "umap_fit", inputs = list(
    space = "pca", dims = 1:3, space_id_out = "umap.pca",
    n_neighbors = 5L, nn_method = "nndescent", seed = 29L))
  promoted <- wb_session_promote_current(b, fit$session, outputs = "model")
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  inputs <- list(model = wb_ref_artifact(model_id, "model"),
    in_space = "pca", dims = 1:3, cell_set = "set:query", min_cells = 2L,
    space_id_out = "umap.pca", overwrite = TRUE)
  first <- wb_session_run_op(wb_session_new(), b, "umap_transform",
    inputs = inputs)
  stale <- wb_session_run_op(wb_session_new(), b, "umap_transform",
    inputs = inputs)
  stopifnot(all(wb_session_promote_current(b, first$session,
    outputs = "projection_membership")$promoted$status == "ok"))
  accepted <- wb_space(b, "umap.pca")
  stale_failure <- try(wb_session_promote_current(b, stale$session,
    outputs = "projection_membership"), silent = TRUE)
  stopifnot(inherits(stale_failure, "try-error"),
    identical(accepted, wb_space(b, "umap.pca")))
  retry <- wb_session_run_op(wb_session_new(), b, "umap_transform",
    inputs = inputs)
  retried <- wb_session_promote_current(b, retry$session,
    outputs = "projection_membership")
  stopifnot(all(retried$promoted$status == "ok"))
  TRUE
}

test_representation_umap_old_core_equivalence_contract <- function() {
  if (!requireNamespace("uwot", quietly = TRUE) ||
      !requireNamespace("rnndescent", quietly = TRUE))
    return(structure(TRUE,
      skip = "uwot and rnndescent are required for the UMAP comparator"))
  path <- if (file.exists(file.path("tests", "fixtures",
      "umap_old_core_reference.json"))) file.path("tests", "fixtures",
        "umap_old_core_reference.json") else file.path("fixtures",
          "umap_old_core_reference.json")
  reference <- jsonlite::read_json(path, simplifyVector = FALSE)
  stopifnot(identical(reference$old_core_commit,
    "651cca3cb0f4b1cd18281cbbee982f267b265e6c"))
  representation_activate_pack()
  root <- tempfile("scwbRepresentation-umap-equivalence-")
  b <- wb_init_bundle(root, overwrite = TRUE)
  on.exit({ try(wb_close_bundle(b), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE) }, add = TRUE)
  cells <- as.character(unlist(reference$cells, use.names = FALSE))
  source <- do.call(rbind, lapply(reference$source, function(x)
    as.numeric(unlist(x, use.names = FALSE))))
  dimnames(source) <- list(cells,
    as.character(unlist(reference$dimensions, use.names = FALSE)))
  wb_define_cells(b, cells); wb_add_space(b, "pca", source)
  packed <- wb_session_run_op(wb_session_new(), b, "umap_fit",
    inputs = list(space = "pca", dims = 1:3, space_id_out = "umap.pca",
      n_neighbors = 5L, min_dist = .2, metric = "euclidean",
      nn_method = "nndescent", spread = 1, seed = 17L))
  expected <- do.call(rbind, lapply(reference$coordinates, function(x)
    as.numeric(unlist(x, use.names = FALSE))))
  dimnames(expected) <- list(cells, c("x", "y"))
  observed <- packed$run$outputs$space_2d$values
  stopifnot(max(abs(observed - expected)) < 1e-12,
    identical(dimnames(observed), dimnames(expected)),
    identical(packed$run$params$dims, 1:3),
    identical(packed$run$params$space_id_out, "umap.pca"),
    identical(packed$run$params$n_neighbors, 5L),
    identical(packed$run$params$min_dist, .2),
    identical(packed$run$params$metric, "euclidean"),
    identical(packed$run$params$nn_method, "nndescent"),
    identical(packed$run$params$spread, 1),
    identical(packed$run$params$seed, 17L))
  TRUE
}

test_representation_umap_contained_model_and_merge_contract <- function() {
  if (!requireNamespace("uwot", quietly = TRUE) ||
      !requireNamespace("rnndescent", quietly = TRUE))
    return(structure(TRUE,
      skip = "uwot and rnndescent are required for the UMAP contract"))
  fixture <- representation_umap_fixture()
  b <- fixture$bundle
  on.exit({ try(wb_close_bundle(b), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE) }, add = TRUE)
  fit <- wb_session_run_op(wb_session_new(), b, "umap_fit", inputs = list(
    space = "pca", dims = 1:3, space_id_out = "umap.pca",
    n_neighbors = 5L, min_dist = .2, metric = "euclidean",
    nn_method = "nndescent", spread = 1, seed = 17L))
  stopifnot(identical(fit$run$op_name, "scwbrepresentation::umap_fit"),
    setequal(names(fit$run$outputs), c("model", "space_2d")),
    identical(dim(fit$run$outputs$space_2d$values), c(30L, 2L)),
    identical(rownames(fit$run$outputs$space_2d$values), fixture$cells),
    identical(fit$run$outputs$model$model_linkage$engine, "uwot"),
    identical(fit$run$outputs$model$model_linkage$model_type, "uwot_umap"))
  promoted <- wb_session_promote_current(b, fit$session, outputs = "model")
  stopifnot(all(promoted$promoted$status == "ok"),
    "umap.pca" %in% wb_spaces(b)$space_id)
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  before <- wb_space(b, "umap.pca")
  wb_close_bundle(b)
  b <- wb_open_bundle(fixture$root)
  fixture$bundle <- b
  model <- wb_art_load(b, model_id, expected_kind = "model_ref")
  stopifnot(file.exists(model$model_path),
    identical(model$model_linkage$engine, "uwot"),
    identical(model$model_linkage$model_type, "uwot_umap"))

  projected <- wb_session_run_op(wb_session_new(), b, "umap_transform",
    inputs = list(model = wb_ref_artifact(model_id, "model"), in_space = "pca", dims = 1:3,
      cell_set = "set:query", min_cells = 2L,
      space_id_out = "umap.pca", overwrite = TRUE))
  stopifnot(identical(projected$run$op_name,
      "scwbrepresentation::umap_transform"),
    identical(rownames(projected$run$outputs$space_2d$values),
      fixture$cells[seq_len(8L)]))
  saved <- wb_session_promote_current(b, projected$session,
    outputs = "projection_membership")
  stopifnot(all(saved$promoted$status == "ok"))
  after <- wb_space(b, "umap.pca")
  before_matrix <- as.matrix(before[, c("x", "y")])
  after_matrix <- as.matrix(after[, c("x", "y")])
  rownames(before_matrix) <- before$cell_id
  rownames(after_matrix) <- after$cell_id
  untouched <- fixture$cells[-seq_len(8L)]
  stopifnot(identical(rownames(after_matrix), fixture$cells),
    identical(after_matrix[untouched, , drop = FALSE],
      before_matrix[untouched, , drop = FALSE]),
    !identical(after_matrix[fixture$cells[seq_len(8L)], , drop = FALSE],
      before_matrix[fixture$cells[seq_len(8L)], , drop = FALSE]),
    isTRUE(wb_bundle_validation_report(fixture$root, "content")$ok))
  TRUE
}

test_representation_umap_corruption_and_atomic_rollback_contract <- function() {
  if (!requireNamespace("uwot", quietly = TRUE) ||
      !requireNamespace("rnndescent", quietly = TRUE))
    return(structure(TRUE,
      skip = "uwot and rnndescent are required for the UMAP contract"))
  fixture <- representation_umap_fixture(); b <- fixture$bundle
  on.exit({ options(scworkbench.space_payload_publication_hook = NULL)
    options(scworkbench.promotion_workspace_publish_hook = NULL)
    try(wb_close_bundle(b), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE) }, add = TRUE)
  fit <- wb_session_run_op(wb_session_new(), b, "umap_fit", inputs = list(
    space = "pca", dims = 1:3, space_id_out = "umap.pca",
    n_neighbors = 5L, nn_method = "nndescent", seed = 19L))
  promoted <- wb_session_promote_current(b, fit$session, outputs = "model")
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  original <- wb_space(b, "umap.pca")
  registry <- DBI::dbGetQuery(b$con,
    "SELECT * FROM space_payloads WHERE space_id='umap.pca';")
  options(scworkbench.space_payload_publication_hook = function(...) {
    stop("injected publication interruption")
  })
  failed <- try(wb_session_run_op(wb_session_new(), b, "umap_transform",
    inputs = list(model = wb_ref_artifact(model_id, "model"), in_space = "pca", dims = 1:3,
      cell_set = "set:query", min_cells = 2L,
      space_id_out = "umap.pca", overwrite = TRUE)), silent = TRUE)
  options(scworkbench.space_payload_publication_hook = NULL)
  stopifnot(inherits(failed, "try-error"),
    identical(original, wb_space(b, "umap.pca")),
    identical(registry, DBI::dbGetQuery(b$con,
      "SELECT * FROM space_payloads WHERE space_id='umap.pca';")))

  retryable <- wb_session_run_op(wb_session_new(), b, "umap_transform",
    inputs = list(model = wb_ref_artifact(model_id, "model"),
      in_space = "pca", dims = 1:3, cell_set = "set:query",
      min_cells = 2L, space_id_out = "umap.pca", overwrite = TRUE))
  before_identity <- scworkbench:::.wb_session_source_identity(b)
  before_payload <- DBI::dbGetQuery(b$con,
    "SELECT * FROM space_payloads WHERE space_id='umap.pca';")
  before_bytes <- digest::digest(file.path(b$dir, before_payload$path[[1L]]),
    algo = "sha256", serialize = FALSE, file = TRUE)
  before_provenance <- DBI::dbGetQuery(b$con,
    "SELECT COUNT(*) AS n FROM prov_runs;")$n[[1L]]
  options(scworkbench.promotion_workspace_publish_hook = function(...) {
    stop("injected durable publication interruption")
  })
  interrupted <- try(wb_session_promote_current(b, retryable$session,
    outputs = "projection_membership"), silent = TRUE)
  options(scworkbench.promotion_workspace_publish_hook = NULL)
  stopifnot(inherits(interrupted, "try-error"),
    identical(before_identity, scworkbench:::.wb_session_source_identity(b)),
    identical(before_payload, DBI::dbGetQuery(b$con,
      "SELECT * FROM space_payloads WHERE space_id='umap.pca';")),
    identical(before_bytes, digest::digest(file.path(b$dir,
      before_payload$path[[1L]]), algo = "sha256", serialize = FALSE,
      file = TRUE)),
    identical(before_provenance, DBI::dbGetQuery(b$con,
      "SELECT COUNT(*) AS n FROM prov_runs;")$n[[1L]]))
  retried <- wb_session_promote_current(b, retryable$session,
    outputs = "projection_membership")
  stopifnot(all(retried$promoted$status == "ok"))
  after_retry <- wb_space(b, "umap.pca")

  model <- wb_art_load(b, model_id, expected_kind = "model_ref")
  path <- model$model_path
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  writeBin(c(bytes, as.raw(0L)), path)
  corrupted <- try(wb_session_run_op(wb_session_new(), b,
    "umap_transform", inputs = list(model = wb_ref_artifact(model_id, "model"), in_space = "pca",
      dims = 1:3, cell_set = "set:query", min_cells = 2L,
      space_id_out = "umap.pca", overwrite = TRUE)), silent = TRUE)
  stopifnot(inherits(corrupted, "try-error"),
    identical(after_retry, wb_space(b, "umap.pca")))
  TRUE
}
