.repr_read_space_axis <- function(reader, capability) {
  chunks <- scworkbench::wb_partition_axis(as.character(capability$cells),
    as.integer(capability$requested_chunk_size))
  blocks <- lapply(chunks, function(cells) reader(cells = cells)$values)
  values <- do.call(rbind, blocks)
  rownames(values) <- as.character(capability$cells)
  values
}

.repr_knn_distances <- function(query, candidate) {
  sqrt(rowSums((candidate - matrix(query, nrow = nrow(candidate),
    ncol = ncol(candidate), byrow = TRUE))^2))
}

.repr_knn_rows <- function(query_id, candidate_id, distance, retain) {
  data.frame(query_cell_id = query_id,
    neighbor_rank = seq_along(retain), cell_id = candidate_id[retain],
    distance = distance[retain], stringsAsFactors = FALSE)
}

.repr_knn_brute_force <- function(query, candidate, k) {
  rows <- lapply(seq_len(nrow(query)), function(i) {
    distance <- .repr_knn_distances(query[i, ], candidate)
    retain <- order(distance, seq_along(distance), method = "radix")
    retain <- retain[seq_len(min(k, length(retain)))]
    .repr_knn_rows(rownames(query)[[i]], rownames(candidate), distance,
      retain)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.repr_knn_rann_exact <- function(query, candidate, k, query_chunk_size) {
  .repr_require("RANN", "exact kNN neighbor search")
  keep_k <- min(as.integer(k), nrow(candidate))
  search_k <- min(nrow(candidate), keep_k + 1L)
  chunks <- scworkbench::wb_partition_axis(rownames(query),
    as.integer(query_chunk_size))
  rows <- list()
  offset <- 0L
  for (cells in chunks) {
    query_index <- match(cells, rownames(query))
    nearest <- RANN::nn2(data = candidate,
      query = query[query_index, , drop = FALSE], k = search_k, eps = 0)
    index_matrix <- nearest$nn.idx
    if (is.null(dim(index_matrix)))
      index_matrix <- matrix(index_matrix, nrow = length(cells), byrow = TRUE)
    index_matrix <- as.matrix(index_matrix)
    if (!identical(as.integer(dim(index_matrix)),
        c(as.integer(length(cells)), as.integer(search_k))))
      stop("RANN returned an invalid exact-neighbor index matrix.",
        call. = FALSE)
    for (j in seq_along(cells)) {
      candidate_index <- as.integer(index_matrix[j, ])
      distance <- .repr_knn_distances(query[query_index[[j]], ],
        candidate[candidate_index, , drop = FALSE])
      ordered <- order(distance, candidate_index, method = "radix")
      candidate_index <- candidate_index[ordered]
      distance <- distance[ordered]
      boundary_tie <- search_k > keep_k &&
        identical(distance[[keep_k]], distance[[keep_k + 1L]])
      if (boundary_tie) {
        full_distance <- .repr_knn_distances(query[query_index[[j]], ],
          candidate)
        retain <- order(full_distance, seq_along(full_distance),
          method = "radix")[seq_len(keep_k)]
        row <- .repr_knn_rows(cells[[j]], rownames(candidate),
          full_distance, retain)
      } else {
        retain <- seq_len(keep_k)
        row <- .repr_knn_rows(cells[[j]], rownames(candidate), distance,
          retain)
        row$cell_id <- rownames(candidate)[candidate_index[retain]]
      }
      offset <- offset + 1L
      rows[[offset]] <- row
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.repr_knn_membership <- function(readers, params, context,
    backend = as.character(params$backend %||% "brute_force")[1L],
    query_chunk_size = as.integer(params$query_chunk_size %||% 1024L)) {
  query <- .repr_read_space_axis(readers$query_space,
    context$typed_space_reads$query_space)
  candidate <- .repr_read_space_axis(readers$candidate_space,
    context$typed_space_reads$candidate_space)
  .repr_assert(ncol(query) == ncol(candidate) && nrow(candidate) > 0L,
    "kNN query and candidate spaces must share dimensions and be nonempty.")
  k <- as.integer(params$k)
  .repr_assert(length(k) == 1L && !is.na(k) && k >= 1L,
    "k must be a positive integer.")
  .repr_assert(backend %in% c("rann_exact", "brute_force"),
    "Unsupported kNN backend.")
  .repr_assert(length(query_chunk_size) == 1L &&
      !is.na(query_chunk_size) && query_chunk_size >= 1L,
    "query_chunk_size must be a positive integer.")
  if (identical(backend, "rann_exact"))
    .repr_knn_rann_exact(query, candidate, k, query_chunk_size) else
    .repr_knn_brute_force(query, candidate, k)
}

.repr_write_membership <- function(writers, rows, context, params) {
  writers$write_terminal_primitive("neighbors",
    scworkbench::wb_prim_obs_table(rows, table_type = "knn_neighbors",
      entity_kind = "cell", meta = list(
        space_id = context$input_ids$space, k = as.integer(params$k),
        include_query = isTRUE(params$include_query))))
  writers$write_terminal_primitive("neighbor_distance",
    scworkbench::wb_prim_cell_vector_family(
      key = rows$query_cell_id, cell_id = rows$cell_id,
      value = rows$distance, value_type = "continuous"))
}

.run_knn_neighbors <- function(readers, writers, params, context) {
  rows <- .repr_knn_membership(readers, params, context,
    backend = as.character(params$backend %||% "rann_exact")[1L],
    query_chunk_size = as.integer(params$query_chunk_size %||% 1024L))
  .repr_write_membership(writers, rows, context, params)
  invisible(list(status = "ok"))
}

.run_knn_pseudobulk <- function(readers, writers, params, context) {
  # Preserve the established pseudobulk scientific implementation. The
  # optimized RANN backend is intentionally scoped to knn_neighbors.
  rows <- .repr_knn_membership(readers, params, context,
    backend = "brute_force")
  features <- as.character(context$typed_slice_reads$matrix_chunks$features)
  candidates <- as.character(
    context$typed_slice_reads$matrix_chunks$cells)
  chunks <- scworkbench::wb_partition_axis(features,
    as.integer(context$typed_slice_reads$matrix_chunks$requested_chunk_size))
  aggregate <- as.character(params$aggregate)[1]
  blocks <- lapply(chunks, function(feature_ids) {
    matrix <- readers$matrix_chunks(features = feature_ids,
      cells = candidates)$values
    do.call(rbind, lapply(as.character(
      context$typed_space_reads$query_space$cells), function(query) {
        neighbor <- rows$cell_id[rows$query_cell_id == query]
        index <- match(neighbor, candidates)
        .repr_assert(!anyNA(index),
          "Accepted kNN membership exceeds the authenticated matrix axis.")
        value <- if (identical(aggregate, "sum"))
          Matrix::rowSums(matrix[, index, drop = FALSE]) else
          Matrix::rowMeans(matrix[, index, drop = FALSE])
        data.frame(feature_id = feature_ids, query_cell_id = query,
          n_neighbors = as.integer(length(index)), aggregate = aggregate,
          value = as.numeric(value), stringsAsFactors = FALSE)
      }))
  })
  # Typed matrix reads are feature-chunked, while the public table contract is
  # query-major. Reassemble deterministically after bounded aggregation.
  table <- do.call(rbind, blocks)
  query_order <- as.character(context$typed_space_reads$query_space$cells)
  table <- table[order(match(table$query_cell_id, query_order),
    match(table$feature_id, features), method = "radix"), , drop = FALSE]
  rownames(table) <- NULL
  writers$write_terminal_primitive("pseudobulk_table",
    scworkbench::wb_prim_feature_table(table,
      table_type = "knn_pseudobulk", meta = list(
        matrix_id = context$input_ids$matrix,
        space_id = context$input_ids$space, k = as.integer(params$k),
        aggregate = aggregate)))
  writers$write_terminal_primitive("neighbor_distance",
    scworkbench::wb_prim_cell_vector_family(
      key = rows$query_cell_id, cell_id = rows$cell_id,
      value = rows$distance, value_type = "continuous"))
  invisible(list(status = "ok"))
}
