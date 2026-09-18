#' Register the provisional scwbRepresentation operation pack
#'
#' Installation is inert. scworkbench invokes this zero-argument entry point
#' only at its explicit operation-pack activation boundary.
#'
#' @return A validated scworkbench operation-pack registration.
#' @export
register_operation_pack <- function() {
  scworkbench::wb_operation_pack_register(
    system.file("scworkbench", "operation-pack.json",
      package = "scwbRepresentation"),
    list(`scwbrepresentation::leiden` = .run_leiden))
}
