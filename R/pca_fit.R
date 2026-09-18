.run_pca_fit <- function(readers, writers, params, context) {
  capability <- context$typed_slice_reads$matrix_chunks
  cells <- as.character(capability$cells)
  features <- as.character(capability$features)
  n_components <- as.integer(params$n_components %||% 50L)
  max_rank <- min(length(cells) - 1L, length(features) - 1L)
  .repr_assert(length(cells) > 2L && length(features) >= 2L,
    "pca_fit requires at least three cells and two selected features")
  .repr_assert(!is.na(n_components) && n_components >= 2L &&
      n_components <= max_rank,
    paste0("pca_fit: n_components=", n_components,
      " exceeds authenticated rank bound ", max_rank, "."))
  chunk_size <- as.integer(params$chunk_size %||% 256L)
  chunks <- scworkbench::wb_partition_axis(features, chunk_size)
  blocks <- lapply(chunks, function(ids) {
    primitive <- readers$matrix_chunks(features = ids, cells = cells)
    t(as.matrix(primitive$values))
  })
  x <- do.call(cbind, blocks)
  dimnames(x) <- list(cells, features)
  approx_gb <- length(cells) * length(features) * 8 / 1024^3
  max_dense_gb <- as.numeric(params$max_dense_gb %||% 8)
  backend <- as.character(params$backend %||% "irlba")[1]
  .repr_assert(is.finite(max_dense_gb) && max_dense_gb > 0,
    "pca_fit: max_dense_gb must be positive")
  .repr_assert(approx_gb <= max_dense_gb || identical(backend, "irlba"),
    "pca_fit: selected block exceeds max_dense_gb")
  feature_mean <- colMeans(x)
  feature_var <- apply(x, 2L, stats::var)
  keep <- is.finite(feature_var) & feature_var > 0
  .repr_assert(sum(keep) >= 2L,
    "pca_fit: fewer than two non-constant features remain")
  fit_features <- features[keep]
  fit_x <- x[, keep, drop = FALSE]
  .repr_assert(n_components <= min(nrow(fit_x) - 1L,
      ncol(fit_x) - 1L),
    "pca_fit: n_components exceeds rank after constant-feature filtering")
  seed <- as.integer(params$seed %||% 1L)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
      inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  fit <- if (identical(backend, "irlba")) {
    .repr_require("irlba", "PCA fitting")
    irlba::prcomp_irlba(fit_x, n = n_components,
      center = isTRUE(params$center), scale. = isTRUE(params$scale))
  } else {
    .repr_assert(identical(backend, "prcomp"),
      "pca_fit: unsupported backend")
    .repr_assert(approx_gb <= max_dense_gb,
      "pca_fit: prcomp exceeds max_dense_gb")
    stats::prcomp(fit_x, rank. = n_components,
      center = isTRUE(params$center), scale. = isTRUE(params$scale))
  }
  labels <- paste0("PC", seq_len(n_components))
  coords <- fit$x[, seq_len(n_components), drop = FALSE]
  dimnames(coords) <- list(cells, labels)
  space_id <- as.character(params$space_id_out)[1]
  writers$write_space("space", structure(list(kind = "space",
    space_id = space_id, values = coords,
    dimensions = as.integer(dim(coords)), cell_id = cells,
    storage_profile = "space_file_backed_hdf5_v1",
    meta = list(source = "pca_fit")), class = "wb_primitive"))
  rotation <- fit$rotation[, seq_len(n_components), drop = FALSE]
  if (is.null(rownames(rotation))) rownames(rotation) <- fit_features
  full_rotation <- matrix(0, nrow = length(features), ncol = n_components,
    dimnames = list(features, labels))
  full_rotation[rownames(rotation), ] <- rotation
  loadings <- data.frame(feature_id = rep(features,
      times = n_components), PC = rep(labels, each = length(features)),
    component = rep(seq_len(n_components), each = length(features)),
    loading = as.numeric(full_rotation), matrix_id = context$input_ids$matrix,
    space_id = space_id, stringsAsFactors = FALSE)
  stdev <- as.numeric(fit$sdev[seq_len(n_components)])
  variance <- stdev^2
  total_variance <- sum(as.numeric(fit$sdev)^2)
  variance_table <- data.frame(PC = labels,
    component = seq_len(n_components), stdev = stdev,
    variance = variance, variance_explained = variance / total_variance,
    cumulative_variance_explained = cumsum(variance / total_variance),
    matrix_id = context$input_ids$matrix, space_id = space_id,
    stringsAsFactors = FALSE)
  summary <- data.frame(matrix_id = context$input_ids$matrix,
    space_id = space_id, n_cells = as.integer(length(cells)),
    n_features_input = as.integer(length(features)),
    n_features = as.integer(length(fit_features)),
    n_features_filtered_constant = as.integer(length(features) -
      length(fit_features)), n_components = n_components,
    center = isTRUE(params$center), scale = isTRUE(params$scale),
    backend = backend, seed = seed, approx_dense_gb = approx_gb,
    max_dense_gb = max_dense_gb, stringsAsFactors = FALSE)
  writers$write_terminal_primitive("loadings",
    scworkbench::wb_prim_feature_table(loadings, table_type = "pca_loadings",
      meta = list(matrix_id = context$input_ids$matrix,
        space_id = space_id)))
  writers$write_terminal_primitive("variance",
    scworkbench::wb_prim_obs_table(variance_table,
      table_type = "pca_variance_explained",
      entity_kind = "principal_component", meta = list(
        matrix_id = context$input_ids$matrix, space_id = space_id)))
  writers$write_terminal_primitive("fit_summary",
    scworkbench::wb_prim_obs_table(summary, table_type = "pca_fit_summary",
      entity_kind = "space", meta = list(
        matrix_id = context$input_ids$matrix, space_id = space_id)))
  invisible(list(status = "ok"))
}
