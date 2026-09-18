.repr_knn_graph <- function(values, k, weight) {
  .repr_require("RANN", "Leiden kNN construction")
  .repr_require("igraph", "Leiden clustering")
  k <- as.integer(k)
  .repr_assert(length(k) == 1L && !is.na(k) && k >= 1L,
    "k must be a positive integer.")
  kk <- min(k + 1L, nrow(values))
  .repr_assert(kk >= 2L, "Not enough cells to compute neighbors.")
  nearest <- RANN::nn2(data = values, query = values, k = kk)
  from <- rep(seq_len(nrow(values)), each = kk - 1L)
  to <- as.vector(t(nearest$nn.idx[, -1L, drop = FALSE]))
  keep <- from != to
  from <- from[keep]
  to <- to[keep]
  a <- pmin(from, to)
  b <- pmax(from, to)
  key <- paste(a, b, sep = "_")
  unique_edge <- !duplicated(key)
  a <- a[unique_edge]
  b <- b[unique_edge]
  graph <- igraph::make_empty_graph(n = nrow(values), directed = FALSE)
  graph <- igraph::add_edges(graph, c(rbind(a, b)))
  igraph::V(graph)$name <- rownames(values)
  if (isTRUE(weight)) {
    distance_from <- rep(seq_len(nrow(values)), each = kk - 1L)
    distance_to <- as.vector(t(nearest$nn.idx[, -1L, drop = FALSE]))
    distance_value <- as.vector(t(nearest$nn.dists[, -1L, drop = FALSE]))
    distance_keep <- distance_from != distance_to
    distance_from <- distance_from[distance_keep]
    distance_to <- distance_to[distance_keep]
    distance_value <- distance_value[distance_keep]
    distance_key <- paste(pmin(distance_from, distance_to),
      pmax(distance_from, distance_to), sep = "_")
    first <- !duplicated(distance_key)
    distance <- data.frame(key = distance_key[first],
      value = distance_value[first], stringsAsFactors = FALSE)
    edge_distance <- distance$value[match(paste(a, b, sep = "_"),
      distance$key)]
    edge_distance[is.na(edge_distance)] <- 1
    igraph::E(graph)$weight <- 1 / (1 + edge_distance)
  }
  graph
}

.run_leiden <- function(readers, writers, params, context) {
  .repr_require("igraph", "Leiden clustering")
  space <- readers$space()
  cells <- as.character(context$effective_cells)
  .repr_assert(length(cells) > 0L && !anyDuplicated(cells),
    "Leiden requires a unique non-empty effective cell axis.")
  values <- .repr_space_matrix(space, cells, params$dims %||% NULL)
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
  graph <- .repr_knn_graph(values, params$k %||% 30L,
    params$weight %||% TRUE)
  use_weights <- isTRUE(params$weight %||% TRUE) &&
    "weight" %in% igraph::edge_attr_names(graph)
  clusters <- igraph::cluster_leiden(graph,
    weights = if (use_weights) igraph::E(graph)$weight else NULL,
    resolution_parameter = as.numeric(params$resolution %||% 1),
    n_iterations = as.integer(params$n_iterations %||% 2L),
    objective_function = "modularity")
  membership <- as.character(as.integer(igraph::membership(clusters)))
  output <- scworkbench::wb_prim_cell_vector(
    stats::setNames(membership, igraph::V(graph)$name),
    value_type = "categorical", meta = list(label = "Leiden clusters"))
  writers$write_terminal_primitive("clusters", output)
  invisible(list(status = "ok"))
}
