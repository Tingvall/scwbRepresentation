args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg),
  ".."), mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
source(file.path(root, "tests", "helper_representation_fixture.R"),
  local = .GlobalEnv)
files <- sort(list.files(file.path(root, "tests"), pattern = "^test_.*\\.R$",
  full.names = TRUE))
for (file in files) source(file, local = .GlobalEnv)
tests <- sort(ls(.GlobalEnv, pattern = "^test_representation_"))
passed <- failed <- skipped <- 0L
for (name in tests) {
  value <- tryCatch(get(name, .GlobalEnv)(), error = identity)
  if (inherits(value, "error")) {
    failed <- failed + 1L
    cat("FAIL", name, conditionMessage(value), "\n")
  } else if (!is.null(attr(value, "skip"))) {
    skipped <- skipped + 1L
    cat("SKIP", name, attr(value, "skip"), "\n")
  } else {
    passed <- passed + 1L
    cat("PASS", name, "\n")
  }
}
cat(sprintf("RESULT discovered=%d executed=%d passed=%d skipped=%d failed=%d\n",
  length(tests), length(tests), passed, skipped, failed))
if (failed) quit(status = 1L)
