test_representation_pca_variance_denominators_contract <- function() {
  representation_activate_pack()
  root <- tempfile("pca-denominators-")
  b <- wb_init_bundle(root)
  on.exit({try(wb_close_bundle(b), silent = TRUE); unlink(root, recursive = TRUE)}, add = TRUE)
  cells <- paste0("c", c(10L, 2L, 1L, setdiff(seq_len(32L), c(10L, 2L, 1L))))
  features <- paste0("g", 1:9)
  set.seed(941)
  x <- matrix(rnorm(length(cells) * 8L), ncol = 8L)
  x <- sweep(x, 2, seq(0.8, 2.2, length.out = 8), `*`)
  x <- sweep(x, 2, seq(1, 3, length.out = 8), `+`)
  x <- cbind(x, constant = 17)
  dimnames(x) <- list(cells, features)
  wb_define_cells(b, cells)
  wb_attach_matrix(b, "input", t(x), feature_type = "gene")
  for (backend in c("prcomp", "irlba")) for (center in c(TRUE, FALSE)) for (scale in c(TRUE, FALSE)) {
    id <- paste("pca", backend, center, scale, sep = "_")
    result <- wb_session_run_op(wb_session_new(), b, "pca_fit", inputs = list(
      matrix = "input", n_components = 2L, backend = backend, center = center,
      scale = scale, seed = 1L, space_id_out = id, chunk_size = 3L))
    v <- result$run$outputs$variance$table
    summary <- result$run$outputs$fit_summary$table
    fitted <- base::scale(x[, 1:8], center = center, scale = scale)
    denominator <- sum(fitted^2) / (nrow(x) - 1L)
    reference <- stats::prcomp(x[, 1:8], center = center, scale. = scale)
    expected <- reference$sdev[1:2]^2
    legacy_total <- if (backend == "prcomp") sum(reference$sdev^2) else sum(expected)
    stopifnot(isTRUE(all.equal(v$variance, expected, tolerance = 1e-5)),
      isTRUE(all.equal(v$variance_explained, expected / legacy_total, tolerance = 1e-5)),
      all(v$variance_explained_denominator == "sum_returned_fit_component_variances"),
      isTRUE(all.equal(v$variance_explained_denominator_value, rep(legacy_total, 2), tolerance = 1e-5)),
      isTRUE(all.equal(v$input_total_variance, rep(denominator, 2), tolerance = 1e-10)),
      isTRUE(all.equal(v$input_total_variance_fraction, expected / denominator, tolerance = 1e-5)),
      abs(tail(v$cumulative_retained_component_variance_fraction, 1) - 1) < 1e-10,
      tail(v$cumulative_input_total_variance_fraction, 1) < 0.99,
      all(v$input_total_variance_status == "available"),
      summary$n_features_input == 9L, summary$n_features == 8L,
      summary$n_features_filtered_constant == 1L,
      abs(summary$input_total_variance - denominator) < 1e-8)
    cat("DENOMINATOR", backend, "center", center, "scale", scale,
      "retained", tail(v$cumulative_retained_component_variance_fraction, 1),
      "input", tail(v$cumulative_input_total_variance_fraction, 1), "\n")
    # Qualify retained public tables for both backends; other settings exercise
    # real fits and actual emitted values without redundant promotions.
    if (center && !scale) {
      promoted <- wb_session_promote_current(b, result$session, outputs = "variance")
      aid <- promoted$promoted$artifact_id[promoted$promoted$output == "variance"]
      stopifnot(all(promoted$promoted$status == "ok"))
      wb_close_bundle(b)
      ro <- wb_open_bundle(root, read_only = TRUE)
      retained <- wb_art_load(ro, aid, expected_kind = "obs_table")$table
      stopifnot(identical(retained, v))
      wb_close_bundle(ro)
      b <- wb_open_bundle(root)
    }
  }
  TRUE
}
