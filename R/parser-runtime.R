.mutt_parser_root <- function() {
  override <- getOption("mutt.parser_root", "")
  if (nzchar(override) && dir.exists(override)) {
    return(normalizePath(override, mustWork = TRUE))
  }

  installed <- system.file("parsers", package = "mutt")
  if (nzchar(installed) && dir.exists(installed)) {
    return(normalizePath(installed, mustWork = TRUE))
  }

  candidates <- c(
    file.path(getwd(), "inst", "parsers"),
    file.path(getwd(), "analysis", "mutt", "inst", "parsers")
  )
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates)) {
    return(normalizePath(candidates[[1L]], mustWork = TRUE))
  }

  stop("MUTT's bundled parser registry could not be found.", call. = FALSE)
}

.mutt_parser_index <- function() {
  root <- .mutt_parser_root()
  index_path <- file.path(root, "index.csv")
  if (file.exists(index_path)) {
    index <- utils::read.csv(index_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!identical(names(index), c("study", "parser_file")) ||
        anyNA(index) || any(!nzchar(index$study)) || any(!nzchar(index$parser_file)) ||
        anyDuplicated(index$study) || anyDuplicated(index$parser_file) ||
        any(!file.exists(file.path(root, index$parser_file)))) {
      stop("MUTT's bundled parser index is invalid.", call. = FALSE)
    }
    return(index)
  }
  ids <- basename(list.dirs(root, recursive = FALSE, full.names = TRUE))
  ids <- ids[file.exists(file.path(root, ids, "parse.R"))]
  data.frame(study = ids, parser_file = file.path(ids, "parse.R"), stringsAsFactors = FALSE)
}

.mutt_parser_ids <- function() {
  .mutt_parser_index()$study
}

.mutt_parser_file <- function(study) {
  index <- .mutt_parser_index()
  row <- match(study, index$study)
  if (is.na(row)) {
    stop("No bundled parser is registered for study `", study, "`.", call. = FALSE)
  }
  path <- file.path(.mutt_parser_root(), index$parser_file[[row]])
  path
}

.mutt_select_studies <- function(studies, available) {
  if (is.list(studies)) studies <- unlist(studies, use.names = TRUE)
  if (is.null(studies)) {
    return(list(studies = available, output_names = available, all_requested = TRUE))
  }
  if (!is.character(studies) || !length(studies) || anyNA(studies) || any(!nzchar(studies))) {
    stop("`studies` must be NULL or a nonempty character vector.", call. = FALSE)
  }
  selected <- unname(studies)
  missing <- setdiff(selected, available)
  if (length(missing)) {
    stop("Unknown study ID(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  output_names <- names(studies)
  if (is.null(output_names)) output_names <- selected
  if (any(!nzchar(output_names)) || anyDuplicated(output_names)) {
    stop("Named `studies` must have unique, nonempty output names.", call. = FALSE)
  }
  list(studies = selected, output_names = output_names, all_requested = FALSE)
}

.mutt_remove_all_na <- function(x) {
  if (is.data.frame(x)) {
    keep <- vapply(x, function(column) !all(is.na(column)), logical(1))
    return(x[, keep, drop = FALSE])
  }
  if (is.list(x)) return(lapply(x, .mutt_remove_all_na))
  x
}

.mutt_run_parser <- function(study, data_root, rawdata, align_samples) {
  parser_environment <- new.env(parent = asNamespace("mutt"))
  sys.source(.mutt_parser_file(study), envir = parser_environment)
  function_name <- paste0("parse_", study)
  if (!exists(function_name, envir = parser_environment, mode = "function", inherits = FALSE)) {
    stop("Bundled parser does not define `", function_name, "()`.", call. = FALSE)
  }

  parser <- get(function_name, envir = parser_environment, inherits = FALSE)
  arguments <- list()
  formal_names <- names(formals(parser))
  if ("raw" %in% formal_names) arguments$raw <- rawdata
  if ("align" %in% formal_names) arguments$align <- align_samples

  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- withCallingHandlers(
    tryCatch(
      suppressMessages(withr::with_dir(data_root, do.call(parser, arguments))),
      error = identity
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  runtime <- proc.time()[["elapsed"]] - started
  list(value = value, warnings = unique(warnings), runtime_seconds = runtime)
}
