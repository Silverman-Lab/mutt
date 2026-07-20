#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: run_functional_study.R STUDY REPOSITORY OUTPUT_DIRECTORY use|rebuild|revalidate",
    call. = FALSE
  )
}

study <- args[[1L]]
repository <- normalizePath(args[[2L]], mustWork = TRUE)
output_directory <- args[[3L]]
mode <- match.arg(tolower(args[[4L]]), c("use", "rebuild", "revalidate"))
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

safe_study <- gsub("[^A-Za-z0-9._-]+", "_", study)
prefix <- file.path(output_directory, safe_study)
started <- Sys.time()
observed_warnings <- character()

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

result <- tryCatch(
  withCallingHandlers(
    mutt::mutt(
      studies = study,
      base_directory = repository,
      align_samples = TRUE,
      functional = switch(
        mode,
        use = TRUE,
        rebuild = "REBUILD",
        revalidate = "REVALIDATE"
      ),
      verbose = TRUE
    ),
    warning = function(condition) {
      observed_warnings <<- c(observed_warnings, conditionMessage(condition))
      message("Warning: ", conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  ),
  error = identity
)

if (inherits(result, "error")) {
  summary <- data.frame(
    study = study,
    parser_status = "error",
    functional_outcome = "error",
    generated = 0L,
    cached = 0L,
    skipped = 0L,
    failed = 0L,
    started_utc = format(started, tz = "UTC", usetz = TRUE),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    warnings = paste(unique(observed_warnings), collapse = " | "),
    error = conditionMessage(result),
    stringsAsFactors = FALSE
  )
  write_tsv(summary, paste0(prefix, "_summary.tsv"))
  quit(save = "no", status = 10L)
}

audit <- attr(result, "audit")
validation <- attr(result, "validation")
provenance <- attr(result, "provenance")
entry <- result[[1L]]
manifest <- if (is.list(entry) && is.list(entry[["function"]])) {
  entry[["function"]]$manifest
} else {
  data.frame()
}

write_tsv(audit, paste0(prefix, "_audit.tsv"))
if (nrow(manifest)) write_tsv(manifest, paste0(prefix, "_functional_manifest.tsv"))

parser_failed <- !identical(audit$status[[1L]], "success")
failed <- if (nrow(manifest)) sum(manifest$status == "failed") else 0L
generated <- if (nrow(manifest)) {
  sum(manifest$status %in% c("generated", "generated_with_warning"))
} else {
  0L
}
cached <- if (nrow(manifest)) sum(manifest$status == "cached") else 0L
skipped <- if (nrow(manifest)) sum(manifest$status == "skipped") else 0L
outcome <- if (parser_failed) {
  "parser_error"
} else if (failed > 0L) {
  "functional_error"
} else if (generated + cached > 0L) {
  "completed"
} else {
  "ineligible"
}

summary <- data.frame(
  study = study,
  parser_status = audit$status[[1L]],
  functional_outcome = outcome,
  generated = generated,
  cached = cached,
  skipped = skipped,
  failed = failed,
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  warnings = paste(unique(c(observed_warnings, audit$warnings)), collapse = " | "),
  error = audit$error[[1L]],
  stringsAsFactors = FALSE
)
write_tsv(summary, paste0(prefix, "_summary.tsv"))
saveRDS(
  list(audit = audit, manifest = manifest, validation = validation, provenance = provenance),
  paste0(prefix, "_diagnostics.rds")
)

print(summary, row.names = FALSE)
if (parser_failed) quit(save = "no", status = 10L)
if (failed > 0L) {
  message("Functional analysis contained ", failed, " failed branch(es).")
  quit(save = "no", status = 11L)
}
