#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(prefix) {
  hit <- args[startsWith(args, paste0("--", prefix, "="))]
  if (length(hit) != 1L) stop("exactly one --", prefix, "= argument is required")
  sub(paste0("--", prefix, "="), "", hit, fixed = TRUE)
}
core <- normalizePath(value("core-source"), mustWork = TRUE)
output <- value("output")
expected_commit <- "19b69bd99bccd81b4e62d839c7d4016984aeec36"
commit <- trimws(system2("git", c("-C", shQuote(core), "rev-parse", "HEAD"),
  stdout = TRUE))
dirty <- system2("git", c("-C", shQuote(core), "status", "--porcelain",
  "--untracked-files=all"), stdout = TRUE)
if (!identical(commit, expected_commit) || any(nzchar(dirty)))
  stop("old-core source must be clean and pinned to ", expected_commit)
devtools::load_all(core, quiet = TRUE)
root <- tempfile("scwbRepresentation-old-core-")
b <- wb_init_bundle(root, overwrite = TRUE)
on.exit({ try(wb_close_bundle(b), silent = TRUE)
  unlink(root, recursive = TRUE, force = TRUE) }, add = TRUE)
cells <- paste0("c", seq_len(8L))
coordinates <- cbind(x = c(0, .1, .2, .1, 10, 10.1, 10.2, 10.1),
  y = c(0, .1, 0, -.1, 10, 10.1, 10, 9.9),
  z = c(1, 1.1, .9, 1, 3, 3.1, 2.9, 3))
rownames(coordinates) <- cells
wb_define_cells(b, cells)
wb_add_space(b, "latent", coordinates, overwrite = TRUE)
result <- wb_session_run_op(wb_session_new(), b, "leiden", inputs = list(
  space = "latent", dims = NULL, k = 5L, resolution = 1,
  n_iterations = 2L, seed = 7L, weight = TRUE, min_cells = 1L))
promotion <- wb_session_promote_current(b, result$session, outputs = "clusters",
  created_by = "phase4a_old_core_reference")
details <- wb_prov_run_details(b, result$run$run_id,
  promotion$provenance_snapshot_id)
identity <- getFromNamespace(".wb_prov_payload_identity", "scworkbench")(
  b, result$run$outputs$clusters)
record <- list(
  fixture_version = "scwbrepresentation-leiden-equivalence-v1",
  old_core_commit = expected_commit,
  input = list(space_id = "latent", effective_cells = cells,
    dimensions = unname(split(coordinates, row(coordinates)))),
  retained_effective_parameters = result$run$params,
  provenance = list(op_name = details$run$op_name,
    symbolic_inputs = details$run$symbolic_inputs,
    execution_inputs = details$run$execution_inputs,
    operation_version = details$run$operation_version),
  terminal = list(kind = result$run$outputs$clusters$kind,
    value_type = result$run$outputs$clusters$value_type,
    cell_id = result$run$outputs$clusters$cell_id,
    value = result$run$outputs$clusters$value,
    meta = result$run$outputs$clusters$meta,
    content_sha256 = identity$content_sha256,
    content_identity_scheme = identity$content_identity_scheme,
    promoted_output_sha256 = as.character(details$outputs$content_sha256[[1L]])))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(record, output, auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = NA)
