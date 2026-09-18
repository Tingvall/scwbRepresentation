#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(prefix) {
  hit <- args[startsWith(args, paste0("--", prefix, "="))]
  if (length(hit) != 1L) stop("exactly one --", prefix,
    "= argument is required")
  sub(paste0("--", prefix, "="), "", hit, fixed = TRUE)
}
core <- normalizePath(value("core-source"), mustWork = TRUE)
output <- value("output")
expected_commit <- "03accf6ac81af0a383cc549d2ae89351af55c7d0"
commit <- trimws(system2("git", c("-C", shQuote(core), "rev-parse", "HEAD"),
  stdout = TRUE))
dirty <- system2("git", c("-C", shQuote(core), "status", "--porcelain",
  "--untracked-files=all"), stdout = TRUE)
if (!identical(commit, expected_commit) || any(nzchar(dirty)))
  stop("old-core source must be clean and pinned to ", expected_commit)
devtools::load_all(core, quiet = TRUE)

root <- tempfile("scwbRepresentation-old-core-pca-")
bundle <- wb_init_bundle(root, overwrite = TRUE)
on.exit({try(wb_close_bundle(bundle), silent = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
cells <- paste0("c", seq_len(12L))
features <- paste0("g", seq_len(5L))
set.seed(17L)
values <- matrix(sample.int(40L, length(cells) * length(features),
  replace = TRUE), nrow = length(features),
  dimnames = list(features, cells))
wb_define_cells(bundle, cells)
wb_attach_matrix(bundle, "RNA", values, feature_type = "gene")
result <- wb_session_run_op(wb_session_new(), bundle, "pca_fit",
  inputs = list(matrix = "RNA", cell_set = "all", feature_ids = features,
    n_components = 3L, center = TRUE, scale = FALSE,
    backend = "prcomp", seed = 1L, space_id_out = "pca",
    overwrite = FALSE, max_dense_gb = 1))
view <- wb_session_view(result$session)
identity <- function(x) {
  out <- scworkbench:::.wb_prov_payload_identity(view, x)
  list(scheme = out$content_identity_scheme, sha256 = out$content_sha256)
}
record <- list(schema = "scwbrepresentation.pca.old_core_reference.v1",
  source_commit = expected_commit, cells = unname(as.list(cells)),
  features = unname(as.list(features)),
  values = unname(lapply(seq_len(nrow(values)), function(i)
    unname(as.list(as.numeric(values[i, ]))))),
  effective_params = result$run$params,
  outputs = lapply(result$run$outputs, identity))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
writeChar(scworkbench:::.wb_operation_pack_canonical_json(record), output,
  eos = NULL, useBytes = TRUE)
