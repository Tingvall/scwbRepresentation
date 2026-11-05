.knn_fixture <- function() {
  representation_activate_pack()
  root <- tempfile("scwbRepresentation-knn-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  cells <- c("q1", "q2", "c1", "c2", "c3")
  features <- c("f1", "f2")
  values <- matrix(c(0, 0, 1, 2, 10, 0, 0, 10, 20, 100),
    nrow = 2L, byrow = TRUE, dimnames = list(features, cells))
  wb_define_cells(bundle, cells)
  wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  coordinates <- cbind(x = c(0, 10, 1, 2, 9), y = rep(0, 5L))
  rownames(coordinates) <- cells
  wb_add_space(bundle, "latent", coordinates, overwrite = TRUE)
  list(root = root, bundle = bundle, cells = cells, features = features)
}

.knn_close <- function(x) {
  try(wb_close_bundle(x$bundle), silent = TRUE)
  unlink(x$root, recursive = TRUE, force = TRUE)
}

.knn_old_core_reference <- function() {
  directory <- if (exists("test_dir", inherits = TRUE))
    get("test_dir", inherits = TRUE) else "tests"
  jsonlite::read_json(file.path(directory, "fixtures",
    "knn_old_core_reference.json"), simplifyVector = TRUE)
}

.knn_content_identity <- function(bundle, value) {
  getFromNamespace(".wb_prov_payload_identity", "scworkbench")(bundle,
    value)
}

test_representation_knn_old_core_equivalence_contract <- function() {
  x <- .knn_fixture(); on.exit(.knn_close(x), add = TRUE)
  reference <- .knn_old_core_reference()
  stopifnot(identical(reference$core_commit,
    "c95f95d5c7fa6b121636408e4cf42f0ab578dd43"))
  neighbors <- wb_session_run_op(wb_session_new(), x$bundle,
    "knn_neighbors", inputs = list(space = "latent",
      query_cells = c("q1", "q2"), candidate_cells = c("c1", "c2", "c3"),
      k = 2L, include_query = FALSE, space_chunk_size = 2L,
      min_cells = 1L))
  expected_table <- data.frame(query_cell_id = rep(c("q1", "q2"), each = 2L),
    neighbor_rank = rep(1:2, 2L), cell_id = c("c1", "c2", "c3", "c2"),
    distance = c(1, 2, 1, 8), stringsAsFactors = FALSE)
  stopifnot(identical(neighbors$run$outputs$neighbors$table, expected_table),
    identical(neighbors$run$outputs$neighbor_distance$key,
      expected_table$query_cell_id),
    identical(neighbors$run$outputs$neighbor_distance$cell_id,
      expected_table$cell_id),
    identical(neighbors$run$outputs$neighbor_distance$value,
      expected_table$distance),
    identical(neighbors$run$params$space_id, "latent"),
    identical(neighbors$run$params$n_query, 2L),
    identical(neighbors$run$params$n_candidates, 3L),
    identical(neighbors$run$params$k, 2L),
    identical(neighbors$run$params$include_query, FALSE),
    identical(.knn_content_identity(x$bundle,
      neighbors$run$outputs$neighbors),
      reference$neighbors$outputs$neighbors$content_identity),
    identical(.knn_content_identity(x$bundle,
      neighbors$run$outputs$neighbor_distance),
      reference$neighbors$outputs$neighbor_distance$content_identity))
  neighbor_promotion <- wb_session_promote_current(x$bundle,
    neighbors$session, outputs = names(neighbors$run$outputs),
    created_by = "phase4e-contract")
  neighbor_details <- wb_prov_run_details(x$bundle, neighbors$run$run_id,
    neighbor_promotion$provenance_snapshot_id)
  stopifnot(all(neighbor_promotion$promoted$status == "ok"),
    all(neighbor_details$outputs$verification_status ==
      "content_verified"),
    all(neighbor_details$outputs$replay_equivalence_policy ==
      "exact_canonical_value"))
  expected <- list(sum = c(3, 30, 12, 120),
    mean = c(1.5, 15, 6, 60))
  for (aggregate in names(expected)) {
    result <- wb_session_run_op(wb_session_new(), x$bundle,
      "knn_pseudobulk", inputs = list(matrix = "RNA", space = "latent",
        query_cells = c("q1", "q2"),
        candidate_cells = c("c1", "c2", "c3"), features = x$features,
        k = 2L, include_query = FALSE, aggregate = aggregate,
        space_chunk_size = 2L, matrix_chunk_size = 1L, min_cells = 1L))
    table <- result$run$outputs$pseudobulk_table$table
    stopifnot(identical(table$feature_id, rep(x$features, 2L)),
      identical(table$query_cell_id, rep(c("q1", "q2"), each = 2L)),
      identical(table$n_neighbors, rep(2L, 4L)),
      identical(table$aggregate, rep(aggregate, 4L)),
      identical(table$value, expected[[aggregate]]),
      identical(result$run$outputs$neighbor_distance$cell_id,
        expected_table$cell_id),
      identical(result$run$params$matrix_id, "RNA"),
      identical(result$run$params$space_id, "latent"),
      identical(result$run$params$n_query, 2L),
      identical(result$run$params$n_candidates, 3L),
      identical(result$run$params$n_features, 2L),
      identical(result$run$params$aggregate, aggregate),
      identical(.knn_content_identity(x$bundle,
        result$run$outputs$pseudobulk_table),
        reference$pseudobulk[[aggregate]]$outputs$pseudobulk_table$
          content_identity),
      identical(.knn_content_identity(x$bundle,
        result$run$outputs$neighbor_distance),
        reference$pseudobulk[[aggregate]]$outputs$neighbor_distance$
          content_identity))
    promoted <- wb_session_promote_current(x$bundle, result$session,
      outputs = names(result$run$outputs), created_by = "phase4e-contract")
    details <- wb_prov_run_details(x$bundle, result$run$run_id,
      promoted$provenance_snapshot_id)
    stopifnot(all(promoted$promoted$status == "ok"),
      all(details$outputs$verification_status == "content_verified"),
      all(details$outputs$replay_equivalence_policy ==
        "exact_canonical_value"))
  }
  TRUE
}

test_representation_knn_overlap_and_bounds_contract <- function() {
  x <- .knn_fixture(); on.exit(.knn_close(x), add = TRUE)
  run <- function(include, k, candidates = x$cells)
    wb_session_run_op(wb_session_new(), x$bundle, "knn_neighbors",
      inputs = list(space = "latent", query_cells = c("q1", "q2"),
        candidate_cells = candidates, k = k, include_query = include,
        min_cells = 1L))$run
  excluded <- run(FALSE, 20L)
  included <- run(TRUE, 20L)
  stopifnot(nrow(excluded$outputs$neighbors$table) == 6L,
    !any(excluded$outputs$neighbors$table$cell_id %in% c("q1", "q2")),
    nrow(included$outputs$neighbors$table) == 10L,
    any(included$outputs$neighbors$table$distance == 0),
    identical(excluded$params$n_candidates, 3L),
    identical(included$params$n_candidates, 5L))
  duplicate <- try(run(TRUE, 2L, c("c1", "c1")), silent = TRUE)
  missing <- try(run(TRUE, 2L, c("c1", "foreign")), silent = TRUE)
  stopifnot(inherits(duplicate, "try-error"), inherits(missing, "try-error"))
  TRUE
}

test_representation_knn_rann_exact_oracle_contract <- function() {
  representation_activate_pack()
  ns <- asNamespace("scwbRepresentation")
  brute <- get(".repr_knn_brute_force", ns)
  rann <- get(".repr_knn_rann_exact", ns)
  candidate <- rbind(c1 = c(0, 0), c2 = c(0, 0), c3 = c(1, 0),
    c4 = c(-1, 0), c5 = c(2, 2), c6 = c(3, 3))
  query <- rbind(q1 = c(0, 0), q2 = c(1, 1), q3 = c(2.5, 2.5),
    q4 = c(-1, 0), q5 = c(10, 10))
  for (k in c(1L, 2L, 4L, 10L)) for (chunk in c(1L, 2L, 5L)) {
    expected <- brute(query, candidate, k)
    observed <- rann(query, candidate, k, chunk)
    stopifnot(identical(observed, expected),
      identical(unique(observed$query_cell_id), rownames(query)))
  }
  tied <- rann(query["q1", , drop = FALSE], candidate, 2L, 1L)
  stopifnot(identical(tied$cell_id, c("c1", "c2")),
    identical(tied$neighbor_rank, 1:2),
    identical(tied$distance, c(0, 0)))
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbRepresentation"),
    simplifyVector = FALSE)
  operation <- manifest$operations[[match("knn_neighbors", vapply(
    manifest$operations, `[[`, character(1), "name"))]]
  dependency <- operation$dependencies[[match("rann_runtime", vapply(
    operation$dependencies, `[[`, character(1), "requirement_id"))]]
  stopifnot(identical(operation$params_schema$properties$backend$default,
      "rann_exact"),
    identical(operation$params_schema$properties$query_chunk_size$default,
      1024L),
    identical(dependency$disposition, "optional"),
    identical(dependency$selected_when$parameter, "backend"),
    identical(dependency$selected_when$equals, "rann_exact"))
  TRUE
}

test_representation_knn_atomic_agreement_contract <- function() {
  x <- .knn_fixture(); on.exit(.knn_close(x), add = TRUE)
  before <- wb_art_list(x$bundle)
  mutate_family <- function(writes, ...) {
    writes$neighbor_distance$value$value[[1]] <-
      writes$neighbor_distance$value$value[[1]] + 1
    writes
  }
  options(scworkbench.operation_pack_v3_output_test_hook = mutate_family)
  failed <- try(wb_session_run_op(wb_session_new(), x$bundle,
    "knn_neighbors", inputs = list(space = "latent", query_cells = "q1",
      candidate_cells = c("c1", "c2"), k = 2L, include_query = FALSE,
      min_cells = 1L)), silent = TRUE)
  options(scworkbench.operation_pack_v3_output_test_hook = NULL)
  stopifnot(inherits(failed, "try-error"), identical(before,
    wb_art_list(x$bundle)))
  omit_table <- function(writes, ...) { writes$pseudobulk_table <- NULL; writes }
  options(scworkbench.operation_pack_v3_output_test_hook = omit_table)
  failed <- try(wb_session_run_op(wb_session_new(), x$bundle,
    "knn_pseudobulk", inputs = list(matrix = "RNA", space = "latent",
      query_cells = "q1", candidate_cells = c("c1", "c2"),
      features = x$features, k = 2L, aggregate = "sum", min_cells = 1L)),
    silent = TRUE)
  options(scworkbench.operation_pack_v3_output_test_hook = NULL)
  stopifnot(inherits(failed, "try-error"), identical(before,
    wb_art_list(x$bundle)))
  TRUE
}

test_representation_knn_pseudobulk_zero_neighbor_contract <- function() {
  x <- .knn_fixture(); on.exit(.knn_close(x), add = TRUE)
  old_hook <- getOption("scworkbench.operation_pack_v3_output_test_hook")
  on.exit(options(scworkbench.operation_pack_v3_output_test_hook =
    old_hook), add = TRUE)
  counts <- function() vapply(c("artifacts", "prov_runs", "prov_commits"),
    function(table) as.integer(DBI::dbGetQuery(x$bundle$con,
      paste0("SELECT COUNT(*) AS n FROM ", table, ";"))$n[[1L]]),
    integer(1))
  before <- list(counts = counts(), artifacts = wb_art_list(x$bundle))
  omit_q2_membership <- function(writes, ...) {
    family <- writes$neighbor_distance$value
    keep <- family$key == "q1"
    family$key <- family$key[keep]
    family$cell_id <- family$cell_id[keep]
    family$value <- family$value[keep]
    writes$neighbor_distance$value <- family
    table <- writes$pseudobulk_table$value$table
    q2 <- table$query_cell_id == "q2"
    table$n_neighbors[q2] <- 0L
    table$value[q2] <- 999
    writes$pseudobulk_table$value$table <- table
    writes
  }
  inputs <- list(matrix = "RNA", space = "latent",
    query_cells = c("q1", "q2"), candidate_cells = c("c1", "c2"),
    features = x$features, k = 1L, include_query = FALSE,
    aggregate = "sum", min_cells = 1L)
  options(scworkbench.operation_pack_v3_output_test_hook =
    omit_q2_membership)
  rejected <- tryCatch(wb_session_run_op(wb_session_new(), x$bundle,
    "knn_pseudobulk", inputs = inputs), error = identity)
  options(scworkbench.operation_pack_v3_output_test_hook = NULL)
  after <- list(counts = counts(), artifacts = wb_art_list(x$bundle))
  stopifnot(inherits(rejected, "wb_operation_pack_error"),
    identical(rejected$code, "PYBENCH_OUTPUT_INVALID"),
    identical(rejected$phase, "audit"),
    grepl("zero_rows_forbidden_v1", conditionMessage(rejected),
      fixed = TRUE), identical(before, after))

  positive <- wb_session_run_op(wb_session_new(), x$bundle,
    "knn_pseudobulk", inputs = inputs)
  table <- positive$run$outputs$pseudobulk_table$table
  stopifnot(all(table$n_neighbors > 0L),
    identical(unique(table$query_cell_id), c("q1", "q2")))
  TRUE
}

test_representation_knn_pack_ownership_contract <- function() {
  representation_activate_pack()
  stopifnot(all(c("scwbrepresentation::knn_neighbors",
    "scwbrepresentation::knn_pseudobulk") %in% wb_list_ops()),
    !exists("wb_op_knn_neighbors", envir = asNamespace("scworkbench"),
      inherits = FALSE),
    !exists("wb_op_knn_pseudobulk", envir = asNamespace("scworkbench"),
      inherits = FALSE))
  TRUE
}
