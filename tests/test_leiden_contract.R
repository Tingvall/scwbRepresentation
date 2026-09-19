test_representation_leiden_scientific_contract <- function() {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("RANN", quietly = TRUE)) return(structure(TRUE,
    skip = "igraph or RANN is not installed"))
  fixture <- representation_fixture()
  on.exit({
    try(wb_close_bundle(fixture$bundle), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  result <- representation_execute(fixture, list(dims = NULL, k = 5L,
    resolution = 1, n_iterations = 2L, seed = 7L, weight = TRUE))
  output <- readRDS(file.path(result$path, "outputs", "clusters",
    "payload.rds"))
  stopifnot(
    identical(output$kind, "cell_vector"),
    identical(output$value_type, "categorical"),
    identical(output$cell_id, fixture$cells),
    identical(output$meta, list(label = "Leiden clusters")),
    length(unique(output$value)) == 2L,
    identical(result$receipt$status, "ok"),
    identical(result$receipt$outputs[[1L]]$name, "clusters"))
  repeat_result <- representation_execute(fixture, list(dims = NULL, k = 5L,
    resolution = 1, n_iterations = 2L, seed = 7L, weight = TRUE))
  repeated <- readRDS(file.path(repeat_result$path, "outputs", "clusters",
    "payload.rds"))
  stopifnot(identical(output, repeated))
  TRUE
}

test_representation_old_core_equivalence_contract <- function() {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("RANN", quietly = TRUE)) return(structure(TRUE,
    skip = "igraph or RANN is not installed"))
  fixture_path <- if (file.exists(file.path("tests", "fixtures",
      "leiden_old_core_reference.json"))) file.path("tests", "fixtures",
      "leiden_old_core_reference.json") else file.path("fixtures",
      "leiden_old_core_reference.json")
  reference <- jsonlite::read_json(fixture_path, simplifyVector = TRUE)
  stopifnot(identical(reference$old_core_commit,
    "19b69bd99bccd81b4e62d839c7d4016984aeec36"))
  fixture <- representation_fixture()
  on.exit({ try(wb_close_bundle(fixture$bundle), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE) }, add = TRUE)
  params <- list(dims = NULL, k = 5L, resolution = 1,
    n_iterations = 2L, seed = 7L, weight = TRUE)
  result <- representation_execute(fixture, params)
  output <- readRDS(file.path(result$path, "outputs", "clusters",
    "payload.rds"))
  identity <- getFromNamespace(".wb_prov_payload_identity", "scworkbench")(
    fixture$bundle, output)
  expected <- list(dims = NULL,
    k = as.integer(reference$retained_effective_parameters$k),
    resolution = as.numeric(reference$retained_effective_parameters$resolution),
    n_iterations = as.integer(reference$retained_effective_parameters$n_iterations),
    seed = as.integer(reference$retained_effective_parameters$seed),
    weight = isTRUE(reference$retained_effective_parameters$weight))
  stopifnot(identical(params, expected),
    identical(as.character(output$cell_id),
      as.character(reference$terminal$cell_id)),
    identical(as.character(output$value),
      as.character(reference$terminal$value)),
    identical(output$meta, reference$terminal$meta),
    identical(identity$content_sha256, reference$terminal$content_sha256),
    identical(identity$content_identity_scheme,
      reference$terminal$content_identity_scheme))
  TRUE
}

test_representation_effective_seed_default_contract <- function() {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("RANN", quietly = TRUE)) return(structure(TRUE,
    skip = "igraph or RANN is not installed"))
  fixture <- representation_fixture()
  on.exit({ try(wb_close_bundle(fixture$bundle), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE) }, add = TRUE)
  common <- list(dims = NULL, k = 5L, resolution = 1,
    n_iterations = 2L, weight = TRUE)
  omitted <- representation_execute(fixture, common)
  explicit <- representation_execute(fixture, c(common, list(seed = 1L)))
  load_output <- function(result) readRDS(file.path(result$path, "outputs",
    "clusters", "payload.rds"))
  set.seed(884L); expected_state <- .Random.seed
  representation_execute(fixture, c(common, list(seed = 19L)))
  observed_after <- stats::runif(3L)
  assign(".Random.seed", expected_state, envir = .GlobalEnv)
  expected_after <- stats::runif(3L)
  stopifnot(identical(load_output(omitted), load_output(explicit)),
    identical(observed_after, expected_after))
  TRUE
}

test_representation_leiden_wrong_input_and_zero_residue_contract <- function() {
  fixture <- representation_fixture()
  on.exit({
    try(wb_close_bundle(fixture$bundle), silent = TRUE)
    unlink(fixture$root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  before <- wb_art_list(fixture$bundle)
  wrong <- wb_prim_cell_vector(stats::setNames(rep("a", 8L), fixture$cells),
    "categorical")
  error <- tryCatch(getFromNamespace(".wb_operation_pack_v2_execute",
    "scworkbench")(fixture$bundle, "scwbrepresentation::leiden",
      list(space = wrong), list(dims = NULL, k = 5L, resolution = 1,
        n_iterations = 2L, seed = 1L, weight = TRUE),
      effective_cells = fixture$cells), error = identity)
  nonfinite <- fixture$space
  nonfinite$values[[1L]] <- Inf
  nonfinite_error <- tryCatch(getFromNamespace(".wb_operation_pack_v2_execute",
    "scworkbench")(fixture$bundle, "scwbrepresentation::leiden",
      list(space = nonfinite), list(dims = NULL, k = 5L, resolution = 1,
        n_iterations = 2L, seed = 1L, weight = TRUE),
      effective_cells = fixture$cells), error = identity)
  misaligned <- fixture$space
  misaligned$values <- misaligned$values[-1L, , drop = FALSE]
  misaligned$cell_id <- misaligned$cell_id[-1L]
  alignment_error <- tryCatch(getFromNamespace(".wb_operation_pack_v2_execute",
    "scworkbench")(fixture$bundle, "scwbrepresentation::leiden",
      list(space = misaligned), list(dims = NULL, k = 5L, resolution = 1,
        n_iterations = 2L, seed = 1L, weight = TRUE),
      effective_cells = fixture$cells), error = identity)
  stopifnot(inherits(error, "error"),
    grepl("Expected kind 'space'", conditionMessage(error), fixed = TRUE),
    inherits(nonfinite_error, "error"),
    grepl("finite", conditionMessage(nonfinite_error), fixed = TRUE),
    inherits(alignment_error, "error"),
    identical(before, wb_art_list(fixture$bundle)))
  TRUE
}

test_representation_manifest_scope_and_dependency_contract <- function() {
  path <- system.file("scworkbench", "operation-pack.json",
    package = "scwbRepresentation")
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  stopifnot(identical(manifest$namespace, "scwbrepresentation"),
    identical(manifest$pack_abi_version, "3"),
    setequal(vapply(manifest$operations, `[[`, "", "name"),
      c("leiden", "pca_fit", "umap_fit", "umap_transform")),
    identical(vapply(manifest$operations[[which(vapply(manifest$operations,
      `[[`, "", "name") == "leiden")]]$dependencies, `[[`, "",
      "name"), c("scworkbench", "igraph", "RANN")))
  TRUE
}
