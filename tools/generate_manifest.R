args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)

identity <- function(name, probe, spec, package = sub("_runtime$", "", name)) list(requirement_id = name,
  kind = "r_package", name = package,
  version = list(scheme = "r_package_version_constraint_v1",
    specifier = spec), disposition = "required", import_probe = probe)
parameter <- function(type, default = NULL, minimum = NULL, maximum = NULL,
                      nullable = is.null(default)) {
  out <- list(type = if (nullable) c(type, "null") else type)
  if (!is.null(default)) out$default <- default
  if (!is.null(minimum)) out$minimum <- minimum
  if (!is.null(maximum)) out$maximum <- maximum
  out
}
params <- list(
  dims = list(type = c("array", "null"), default = NULL,
    uniqueItems = TRUE, items = list(type = "integer", minimum = 1L)),
  k = parameter("integer", 30L, 5L, 200L, FALSE),
  resolution = parameter("number", 1, .01, 10, FALSE),
  n_iterations = parameter("integer", 2L, 1L, nullable = FALSE),
  seed = parameter("integer", 1L, nullable = FALSE),
  weight = list(type = "boolean", default = TRUE))
retained <- lapply(names(params), function(name) {
  out <- list(name = name, source = "parameter_value", parameter = name)
  if (identical(name, "dims")) out$required <- FALSE
  out
})
manifest <- list(
  schema_version = "scworkbench.operation_pack.v3", pack_abi_version = "3",
  identity_scheme = "scworkbench.operation_pack.canonical_json.v1",
  namespace = "scwbrepresentation",
  package = list(name = "scwbRepresentation", version = "0.1.0"),
  backend = "R", compatibility = list(
    scworkbench = list(minimum = "0.1.0", maximum_exclusive = "0.2.0"),
    pack_abi = list(minimum = "3", maximum = "3")),
  entry_point = "scwbRepresentation::register_operation_pack",
  operations = list(list(
    name = "leiden", qualified_name = "scwbrepresentation::leiden",
    title = "Leiden clustering",
    description = "Graph clustering on a kNN graph derived from an embedding space.",
    runner = ".run_leiden",
    inputs = list(list(name = "space", kind = "space", required = TRUE,
      alignment = "cell_ordered")),
    params_schema = list(type = "object", properties = params,
      additionalProperties = FALSE),
    dimension_policies = list(list(output = "clusters",
      policy = "authenticated_output_axes_v1", shape = "one_per_cell",
      cell_axis = list(source = "effective_cells_v1", input = "space"),
      row_order = "input_cell_order_v1", missing_values = "forbidden")),
    outputs = list(list(name = "clusters", category = "terminal_artifact",
      kind = "cell_vector", storage_profile = "terminal_primitive_v1",
      allocated_id_key = "clusters_id", required = TRUE,
      promotion_group = "clusters", value_kind = "categorical")),
    dependencies = list(
      identity("scworkbench_runtime", "wb_prim_cell_vector", ">=0.1.0,<0.2.0"),
      identity("igraph_runtime", "cluster_leiden", ">=2.0.0,<3.0.0"),
      identity("rann_runtime", "nn2", ">=2.6.0,<3.0.0", "RANN")),
    promotion_closure = list(clusters = list("clusters")),
    replay_equivalence = list(clusters = "exact_canonical_value"),
    retained_parameters = retained,
    public_interface = list(intent = list("cluster"),
      inputs = list(list(name = "space", label = "space",
        help = "Cell-aligned embedding space.")),
      parameters = lapply(names(params), function(name) list(name = name,
        type = if (name == "dims") "list" else if (name %in%
          c("k", "n_iterations", "seed")) "int" else if (name == "weight")
          "bool" else "dbl", help = paste("Leiden", name, "parameter."))),
      outputs = list(list(name = "clusters",
        description = "Categorical Leiden label for each effective cell."))))))
path <- file.path(root, "inst", "scworkbench", "operation-pack.json")
jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = NA)
