`%||%` <- function(x, y) if (is.null(x)) y else x

.repr_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.repr_require <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s.", package, purpose),
      call. = FALSE)
  }
}

.repr_space_matrix <- function(space, cells, dims = NULL) {
  .repr_assert(inherits(space, "wb_primitive") &&
      identical(as.character(space$kind)[1L], "space"),
    "Expected a space primitive.")
  values <- space$values
  .repr_assert(is.matrix(values) && !is.null(rownames(values)),
    "space must carry a coordinate matrix with cell row names.")
  cells <- as.character(cells)
  index <- match(cells, rownames(values))
  .repr_assert(!anyNA(index),
    "Effective cells are absent from the supplied space primitive.")
  values <- values[index, , drop = FALSE]
  rownames(values) <- cells
  if (!is.null(dims)) {
    dims <- as.integer(unlist(dims, use.names = FALSE))
    .repr_assert(length(dims) > 0L && !anyDuplicated(dims) &&
        all(dims >= 1L & dims <= ncol(values)),
      "dims must be unique in-range one-based space dimensions.")
    values <- values[, dims, drop = FALSE]
  }
  .repr_assert(ncol(values) >= 1L && all(is.finite(values)),
    "Space coordinates must be finite and non-empty.")
  values
}
