#!/usr/bin/env Rscript

library(mutt)

data_root <- Sys.getenv("MUTT_DATA_DIR", unset = "")
output_dir <- Sys.getenv("MUTT_CHECK_OUTPUT", unset = "")
align_value <- tolower(trimws(Sys.getenv("MUTT_CHECK_ALIGN_SAMPLES", unset = "false")))
if (!align_value %in% c("false", "true")) {
  stop("MUTT_CHECK_ALIGN_SAMPLES must be 'true' or 'false'.", call. = FALSE)
}
align_samples <- identical(align_value, "true")
if (!nzchar(data_root) || !dir.exists(data_root)) {
  stop("MUTT_DATA_DIR must name the studies directory.", call. = FALSE)
}
if (!nzchar(output_dir)) stop("MUTT_CHECK_OUTPUT is required.", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

index_path <- system.file("parsers", "index.csv", package = "mutt")
index <- utils::read.csv(index_path, stringsAsFactors = FALSE, check.names = FALSE)
study_ids <- index$study
report_path <- file.path(output_dir, "all_parser_check.csv")
report <- if (file.exists(report_path)) {
  utils::read.csv(report_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame()
}
if (nrow(report) && !"align_samples" %in% names(report)) report$align_samples <- FALSE
if (nrow(report) && any(report$align_samples != align_samples)) {
  stop("The existing report was created with a different alignment mode.", call. = FALSE)
}
requested <- Sys.getenv("MUTT_CHECK_STUDIES", unset = "")
if (nzchar(requested)) {
  requested <- trimws(strsplit(requested, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(requested, study_ids)
  if (length(unknown)) stop("Unknown requested studies: ", paste(unknown, collapse = ", "), call. = FALSE)
  study_ids <- requested
  if (nrow(report)) report <- report[!report$study %in% study_ids, , drop = FALSE]
}
completed <- if (nrow(report)) report$study else character()

for (i in seq_along(study_ids)) {
  study <- study_ids[[i]]
  if (study %in% completed) next
  cat(sprintf("[%d/%d] %s\n", i, length(study_ids), study))
  emitted <- character()
  started <- Sys.time()
  result <- withCallingHandlers(
    tryCatch(
      mutt(
        studies = study,
        base_directory = data_root,
        align_samples = align_samples,
        functional = FALSE
      ),
      error = identity
    ),
    warning = function(w) {
      emitted <<- c(emitted, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(result, "error")) {
    row <- data.frame(
      study = study,
      status = "error",
      align_samples = align_samples,
      runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      data_source = "",
      nonempty_entries = 0L,
      warnings = paste(unique(emitted), collapse = " | "),
      error = conditionMessage(result),
      stringsAsFactors = FALSE
    )
  } else {
    audit <- attr(result, "audit")
    object <- result[[1L]]
    nonempty <- if (is.list(object)) {
      sum(!vapply(object, function(x) is.null(x) || (length(x) == 1L && all(is.na(x))), logical(1)))
    } else {
      0L
    }
    row <- data.frame(
      study = study,
      status = audit$status[[1L]],
      align_samples = align_samples,
      runtime_seconds = audit$runtime_seconds[[1L]],
      data_source = audit$data_source[[1L]],
      nonempty_entries = nonempty,
      warnings = paste(unique(c(emitted, audit$warnings[[1L]])), collapse = " | "),
      error = audit$error[[1L]],
      stringsAsFactors = FALSE
    )
  }

  report <- rbind(report, row)
  utils::write.csv(report, report_path, row.names = FALSE, na = "")
  cat(sprintf("  %s (%.1f s)\n", row$status, row$runtime_seconds))
  rm(result)
  invisible(gc())
}

utils::write.csv(report, report_path, row.names = FALSE, na = "")
summary <- aggregate(study ~ align_samples + status, report, length)
names(summary)[[3L]] <- "n"
utils::write.csv(summary, file.path(output_dir, "all_parser_check_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
print(summary, row.names = FALSE)
if (any(report$status != "success")) quit(save = "no", status = 1L)
