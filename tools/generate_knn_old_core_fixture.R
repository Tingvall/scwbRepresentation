#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
prefix <- "--core-source="
source_arg <- args[startsWith(args, prefix)]
if (length(source_arg) != 1L) stop("one --core-source argument is required")
core <- normalizePath(sub(prefix, "", source_arg, fixed = TRUE),
  winslash = "/", mustWork = TRUE)
expected <- "c95f95d5c7fa6b121636408e4cf42f0ab578dd43"
commit <- system2("git", c("-C", shQuote(core), "rev-parse", "HEAD"),
  stdout = TRUE)
if (!identical(trimws(commit), expected)) stop("old-core commit mismatch")
dirty <- system2("git", c("-C", shQuote(core), "status", "--porcelain"),
  stdout = TRUE)
if (length(dirty) && any(nzchar(dirty))) stop("old-core source is dirty")
devtools::load_all(core, quiet = TRUE)
root <- tempfile("phase4e-old-core-")
bundle <- wb_init_bundle(root, overwrite = TRUE)
on.exit({try(wb_close_bundle(bundle), silent = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
cells <- c("q1", "q2", "c1", "c2", "c3")
features <- c("f1", "f2")
values <- matrix(c(0, 0, 1, 2, 10, 0, 0, 10, 20, 100), nrow = 2L,
  byrow = TRUE, dimnames = list(features, cells))
wb_define_cells(bundle, cells)
wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
  feature_type = "gene")
coordinates <- cbind(x = c(0, 10, 1, 2, 9), y = rep(0, 5L))
rownames(coordinates) <- cells
wb_add_space(bundle, "latent", coordinates, overwrite = TRUE)
arguments <- list(space_id = "latent", query_cells = c("q1", "q2"),
  candidate_cells = c("c1", "c2", "c3"), k = 2L,
  include_query = FALSE)
neighbors <- do.call(wb_op_knn_neighbors, c(list(bundle = bundle), arguments))
pseudobulk <- lapply(c("sum", "mean"), function(aggregate)
  do.call(wb_op_knn_pseudobulk, c(list(bundle = bundle, matrix_id = "RNA",
    features = features, aggregate = aggregate), arguments)))
names(pseudobulk) <- c("sum", "mean")
primitive <- function(value) list(value = unclass(value),
  content_identity = scworkbench:::.wb_prov_payload_identity(bundle, value))
record <- list(
  schema = "scwbrepresentation.knn.old_core_reference.v1",
  core_commit = expected,
  inputs = list(cells = cells, features = features,
    values = unname(lapply(seq_len(nrow(values)), function(i)
      unname(as.list(values[i, ])))), coordinates = unname(lapply(
        seq_len(nrow(coordinates)), function(i)
          unname(as.list(coordinates[i, ]))))),
  arguments = arguments,
  neighbors = list(params = neighbors$params,
    outputs = lapply(neighbors$outputs, primitive)),
  pseudobulk = lapply(pseudobulk, function(result) list(params = result$params,
    outputs = lapply(result$outputs, primitive))))
script <- sub("^--file=", "", commandArgs(FALSE)[grepl(
  "^--file=", commandArgs(FALSE))][1])
pack_root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
path <- file.path(pack_root, "tests", "fixtures",
  "knn_old_core_reference.json")
jsonlite::write_json(record, path, auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = NA)
cat("WROTE", path, "\n")
