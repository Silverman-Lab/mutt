#' Parse MUTT microbiome studies
#'
#' `mutt()` is the package's public API. It resolves study data from a local
#' Git LFS checkout or the versioned MUTT user cache, runs the bundled study
#' parsers in isolation, validates their returned objects, and optionally runs
#' PICRUSt2 for eligible ASV count tables and FAPROTAX for eligible taxonomic
#' abundance tables. FAPROTAX prefers counts and falls back to proportions only
#' when a matching count table is unavailable.
#'
#' @param studies `NULL`, a character vector of study IDs, or a named character
#'   vector. `NULL` parses every locally available registered study. Names, when
#'   supplied, become names in the returned list.
#' @param base_directory Optional local data root. This may be the MUTT
#'   repository root containing `studies/` or the `studies/` directory itself.
#'   When omitted, MUTT checks `MUTT_DATA_DIR`, `./studies`, and its versioned
#'   user cache.
#' @param rawdata Logical scalar passed to parsers as `raw`.
#' @param align_samples Logical scalar passed to parsers as `align`.
#' @param save_to Optional path ending in `.RData`. The parsed result and a
#'   companion validation object are saved without changing the returned value.
#' @param verbose Logical scalar controlling progress and parser summaries.
#' @param functional `FALSE`, `TRUE`, or `"REBUILD"`. `TRUE` runs or reuses
#'   PICRUSt2 and FAPROTAX where valid inputs and tools are available;
#'   `"REBUILD"` forces eligible cached branches to rerun.
#'
#' @return A `mutt_result`: a named list whose successful entries retain the
#'   established MUTT study structure. `attr(x, "audit")` contains one row per
#'   requested study, `attr(x, "provenance")` records package and data-release
#'   information, and `attr(x, "validation")` contains structural checks.
#'
#' @examples
#' \dontrun{
#' one <- mutt(
#'   studies = "2021_vandeputte_naturecommunications_flow_timeseries",
#'   base_directory = "/path/to/mutt"
#' )
#'
#' with_functions <- mutt(
#'   studies = "2022_cvandevelde_ismecommunications_culturedflowhumanfecal",
#'   base_directory = "/path/to/mutt/studies",
#'   functional = TRUE
#' )
#' }
#'
#' @export
mutt <- function(
  studies = NULL,
  base_directory = NULL,
  rawdata = FALSE,
  align_samples = FALSE,
  save_to = NULL,
  verbose = FALSE,
  functional = FALSE
) {
  .mutt_assert_flag(rawdata, "rawdata")
  .mutt_assert_flag(align_samples, "align_samples")
  .mutt_assert_flag(verbose, "verbose")
  if (!is.null(save_to) &&
      (!is.character(save_to) || length(save_to) != 1L || is.na(save_to) || !nzchar(save_to))) {
    stop("`save_to` must be NULL or one nonempty path.", call. = FALSE)
  }

  functional_mode <- .functional_mode(functional)
  parser_ids <- .mutt_parser_ids()
  selection <- .mutt_select_studies(studies, parser_ids)
  selected <- selection$studies
  output_names <- selection$output_names
  locations <- .mutt_resolve_study_locations(
    selected,
    base_directory = base_directory,
    all_requested = selection$all_requested
  )
  if (rawdata) {
    missing_references <- any(vapply(
      locations,
      function(location) !dir.exists(file.path(location$data_root, "helperdata")),
      logical(1)
    ))
    if (missing_references) {
      stop(
        "`rawdata = TRUE` requires the reference files from a local Git LFS checkout. ",
        "Set `base_directory` or MUTT_DATA_DIR to that checkout.",
        call. = FALSE
      )
    }
  }

  functional_tools <- if (functional_mode == "off") NULL else .discover_functional_tools()
  if (functional_mode != "off") {
    if (!nzchar(functional_tools$picrust2)) {
      message("PICRUSt2 was not found on PATH; eligible branches will be recorded as skipped.")
    }
    if (!nzchar(functional_tools$faprotax) || !nzchar(functional_tools$faprotax_db)) {
      message("FAPROTAX was not found with its database; eligible branches will be recorded as skipped.")
    }
  }

  result <- setNames(vector("list", length(selected)), output_names)
  validation <- setNames(vector("list", length(selected)), output_names)
  audit_rows <- vector("list", length(selected))
  progress <- NULL
  if (!verbose && length(selected) > 1L) {
    progress <- utils::txtProgressBar(min = 0, max = length(selected), style = 3)
    on.exit(try(close(progress), silent = TRUE), add = TRUE)
  }

  for (index in seq_along(selected)) {
    study <- selected[[index]]
    output_name <- output_names[[index]]
    location <- locations[[study]]
    if (verbose) message(sprintf("[%d/%d] %s", index, length(selected), study))

    parsed <- .mutt_run_parser(
      study = study,
      data_root = location$data_root,
      rawdata = rawdata,
      align_samples = align_samples
    )
    status <- "success"
    error_message <- ""
    parser_warnings <- parsed$warnings

    if (inherits(parsed$value, "error")) {
      status <- "error"
      error_message <- conditionMessage(parsed$value)
    } else if (!is.list(parsed$value)) {
      status <- "error"
      error_message <- "Parser did not return a list."
    } else {
      value <- .mutt_remove_all_na(parsed$value)
      if (!is.null(value$tax)) value$tax <- add_sequence_column(value$tax)
      functional_value <- tryCatch(
        .run_study_functional(
          parsed = value,
          study_dir = location$study_dir,
          mode = functional_mode,
          tools = functional_tools
        ),
        error = identity
      )
      if (inherits(functional_value, "error")) {
        parser_warnings <- c(
          parser_warnings,
          paste("Functional inference failed:", conditionMessage(functional_value))
        )
        functional_value <- .empty_functional_result(functional_mode != "off")
      }
      value[["function"]] <- functional_value
      value <- standardize_output_order(value)
      if (exists("datasets", inherits = TRUE) && !is.null(datasets[[study]])) {
        value$studydemographics <- datasets[[study]]
      }
      result[[output_name]] <- value
      validation[[output_name]] <- validate_output_structure(value, study_name = study)
    }

    audit_rows[[index]] <- data.frame(
      output_name = output_name,
      study = study,
      status = status,
      data_source = location$source,
      data_directory = location$study_dir,
      runtime_seconds = unname(parsed$runtime_seconds),
      warnings = paste(unique(parser_warnings), collapse = " | "),
      error = error_message,
      stringsAsFactors = FALSE
    )
    if (!is.null(progress)) utils::setTxtProgressBar(progress, index)
  }
  if (!is.null(progress)) close(progress)

  audit <- do.call(rbind, audit_rows)
  provenance <- list(
    package_version = tryCatch(as.character(utils::packageVersion("mutt")), error = function(e) "development"),
    data_release = .mutt_data_release,
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    functional_mode = functional_mode
  )
  result <- .new_mutt_result(result, audit, provenance, validation)

  if (!is.null(save_to)) .mutt_save_result(result, validation, save_to, verbose)
  failed <- audit$status != "success"
  if (any(failed)) {
    warning(
      sum(failed), " of ", nrow(audit),
      " requested study parsers failed. Inspect attr(result, \"audit\").",
      call. = FALSE
    )
  }
  result
}

.mutt_assert_flag <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(value)
}

.mutt_save_result <- function(result, validation, save_to, verbose) {
  directory <- dirname(save_to)
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  microbialscalerepository <- result
  save(microbialscalerepository, file = save_to)
  validation_file <- if (grepl("\\.RData$", save_to, ignore.case = TRUE)) {
    sub("\\.RData$", "_validation.RData", save_to, ignore.case = TRUE)
  } else {
    paste0(save_to, "_validation.RData")
  }
  validation_results <- validation
  save(validation_results, file = validation_file)
  if (verbose) {
    message("Saved result to ", normalizePath(save_to, mustWork = TRUE))
    message("Saved validation to ", normalizePath(validation_file, mustWork = TRUE))
  }
  invisible(save_to)
}
