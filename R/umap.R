.repr_umap_preflight <- function(method) {
  method <- as.character(method %||% "nndescent")[1]
  .repr_assert(method %in% c("nndescent", "annoy", "fnn"),
    "umap_fit: unsupported nn_method")
  if (identical(method, "nndescent"))
    .repr_require("rnndescent", "UMAP nndescent neighbor search")
  if (identical(method, "fnn"))
    .repr_require("FNN", "UMAP FNN neighbor search")
  method
}

.repr_selected_space <- function(space, cells, dims = NULL) {
  values <- as.matrix(space$values)
  index <- match(cells, rownames(values))
  .repr_assert(!anyNA(index),
    "UMAP source space does not cover every effective cell")
  values <- values[index, , drop = FALSE]
  dims <- as.integer(unlist(dims %||% seq_len(ncol(values)),
    use.names = FALSE))
  .repr_assert(length(dims) && !anyNA(dims) && all(dims >= 1L) &&
      all(dims <= ncol(values)) && !anyDuplicated(dims),
    "UMAP dimensions are empty, duplicated, or out of bounds")
  values[, dims, drop = FALSE]
}

.run_umap_fit <- function(readers, writers, params, context) {
  .repr_require("uwot", "UMAP fitting")
  space <- readers$space()
  cells <- as.character(context$effective_cells)
  x <- .repr_selected_space(space, cells, params$dims %||% NULL)
  method <- .repr_umap_preflight(params$nn_method)
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
  fit <- uwot::umap(x,
    n_neighbors = as.integer(params$n_neighbors %||% 30L),
    min_dist = as.numeric(params$min_dist %||% .3),
    metric = as.character(params$metric %||% "euclidean")[1],
    nn_method = method, spread = as.numeric(params$spread %||% 1),
    n_threads = 1L, n_sgd_threads = 1L, seed = seed,
    ret_model = TRUE, verbose = FALSE)
  coordinates <- fit$embedding
  .repr_assert(is.matrix(coordinates) && ncol(coordinates) == 2L &&
      all(is.finite(coordinates)), "uwot returned invalid coordinates")
  rownames(coordinates) <- cells
  colnames(coordinates) <- c("x", "y")
  target <- as.character(params$space_id_out)[1]
  writers$write_space_2d("space_2d", structure(list(kind = "space_2d",
    space_id = target, values = coordinates,
    dimensions = as.integer(dim(coordinates)), cell_id = cells,
    storage_profile = "space_2d_dense_v1",
    meta = list(source = "uwot_fit")), class = "wb_primitive"))
  tree <- writers$open_model_tree("model")
  uwot::save_uwot(fit, tree$path(tree$primary_path), verbose = FALSE)
  model <- scworkbench::wb_prim_model_ref(target,
    as.character(context$input_ids$space)[1], engine = "uwot")
  model$model_type <- "uwot_umap"
  tree$finish(model)
  invisible(list(status = "ok"))
}

.run_umap_transform <- function(readers, writers, params, context) {
  .repr_require("uwot", "UMAP transformation")
  payload <- readers$fitted_model()
  model <- uwot::load_uwot(payload$path(payload$primary_path),
    verbose = FALSE)
  space <- readers$in_space()
  cells <- as.character(context$effective_cells)
  x <- .repr_selected_space(space, cells, params$dims %||% NULL)
  coordinates <- uwot::umap_transform(x, model, verbose = FALSE)
  .repr_assert(is.matrix(coordinates) && ncol(coordinates) == 2L &&
      all(is.finite(coordinates)),
    "uwot transform returned invalid coordinates")
  rownames(coordinates) <- cells
  colnames(coordinates) <- c("x", "y")
  target <- as.character(params$space_id_out)[1]
  writers$write_space_2d("space_2d", structure(list(kind = "space_2d",
    space_id = target, values = coordinates,
    dimensions = as.integer(dim(coordinates)), cell_id = cells,
    storage_profile = "space_2d_dense_v1",
    meta = list(source = "uwot_transform")), class = "wb_primitive"))
  writers$write_terminal_primitive("projection_membership",
    scworkbench::wb_prim_cell_vector(stats::setNames(
      rep("projected", length(cells)), cells), value_type = "categorical",
      meta = list(space_id = target, model_bound = TRUE)))
  invisible(list(status = "ok"))
}
