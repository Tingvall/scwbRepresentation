.repr_read_space_axis <- function(reader, capability) {
  chunks <- scworkbench::wb_partition_axis(as.character(capability$cells),
    as.integer(capability$requested_chunk_size))
  blocks <- lapply(chunks, function(cells) reader(cells = cells)$values)
  values <- do.call(rbind, blocks)
  rownames(values) <- as.character(capability$cells)
  values
}

.repr_knn_membership <- function(readers, params, context) {
  query <- .repr_read_space_axis(readers$query_space,
    context$typed_space_reads$query_space)
  candidate <- .repr_read_space_axis(readers$candidate_space,
    context$typed_space_reads$candidate_space)
  .repr_assert(ncol(query) == ncol(candidate) && nrow(candidate) > 0L,
    "kNN query and candidate spaces must share dimensions and be nonempty.")
  k <- as.integer(params$k)
  rows <- vector("list", nrow(query))
  for (i in seq_len(nrow(query))) {
    distance <- sqrt(rowSums((candidate - matrix(query[i, ],
      nrow = nrow(candidate), ncol = ncol(candidate), byrow = TRUE))^2))
    order_index <- order(distance, seq_along(distance), method = "radix")
    retain <- order_index[seq_len(min(k, length(order_index)))]
    rows[[i]] <- data.frame(query_cell_id = rownames(query)[[i]],
      neighbor_rank = seq_along(retain),
      cell_id = rownames(candidate)[retain], distance = distance[retain],
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
  rows <- .repr_knn_membership(readers, params, context)
  .repr_write_membership(writers, rows, context, params)
  invisible(list(status = "ok"))
}

.run_knn_pseudobulk <- function(readers, writers, params, context) {
  rows <- .repr_knn_membership(readers, params, context)
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
