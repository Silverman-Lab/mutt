# ─────────────────────────────────────────────────────────────────────────────
#' Microbial Scale Repository
#'
#' Public genomic (Amplicon or Shotgun Sequencing) datasets parsed or reanalyzed
#' stored in subdirectories under `base_directory` that contain a `parse.R`
#' file, source each in isolation, run its `parse_<folder>()` function, and
#' collect the results. 
#'
#' @author Maxwell Konnaris
#' @author Manan Saxena
#' @author Justin Silverman
#' @maintainer Justin Silverman <justinsilverman@psu.edu>
#' @maintainer Maxwell Konnaris <mak6930@psu.edu>
#'
#' @importFrom utils txtProgressBar setTxtProgressBar 
#' @import tidyverse 
#' @import optparse
#'
#' @param studies Character vector *or named list* indicating which parser
#'        subdirectories to run.  If unnamed, `studies` entries must exactly
#'        match folder names.  If named, `names(studies)` become the output
#'        list names and the values are folder names.  `NULL` (default) runs
#'        every parser folder found.
#' @param base_directory Character(1). Directory containing the parser
#'        sub‑folders.  Default `"."`.
#' @param rawdata Logical(1). Passed as `raw = rawdata` to each parser.
#'        Default `FALSE`.
#' @param align_samples Logical(1). If `TRUE`, aligns the rows of `$counts`
#'        to the sample IDs in the first column of `$scale`.  Default
#'        `FALSE`.
#' @param save_to Optional character(1) path to save the resulting list as an
#'        `.RData` file.
#'
#' @return Named list of parser outputs.  Each element/study must at minimum contain
#'         a `$counts` data.frame and a `$scale` data.frame. Details listed below
#'         in section.
#' 
#' @section Individual Parser Output Specification:
#'
#' Each `parse.R` must define a single function named `parse_<foldername>` where
#' `<foldername>` is the lowercase name of the parser subdirectory. The function
#' should return a named list with the elements below. Each can be a list representing multiple 
#' datasets (e.g. both amplicon and shotgun data) either from the original manuscript or reanalyzed data.
#'
#' \describe{
#'   \item{\code{counts}}{
#'     An integer-valued matrix of raw counts (not a data.frame), dimension \eqn{N \times D},
#'     where rows are sample IDs and columns are sequence IDs (e.g., taxa IDs).
#'     Must include sample ID rownames that link to \code{proportions}, \code{scale},
#'     and \code{metadata}.
#'   }
#'
#'   \item{\code{proportions}}{
#'     A real-valued matrix of relative abundances (not a data.frame),
#'     dimension \eqn{N \times D}, with row and column names matching sample and
#'     sequence IDs respectively. Must include sample ID rownames linking to
#'     \code{counts}, \code{scale}, and \code{metadata}.
#'   }
#'
#'   \item{\code{scale}}{
#'     A matrix of positive-valued scaling factors (often \eqn{N \times 1}, but
#'     may vary depending on method, e.g., with columns for mean and SD).
#'     Must include sample ID rownames matching other tables.
#'   }
#'
#'   \item{\code{metadata} (optional)}{
#'     A data.frame of dimension \eqn{N \times Q} containing metadata for each sample.
#'     Must include a column or rownames matching sample IDs in \code{counts} and \code{proportions}.
#'   }
#'
#'   \item{\code{tax} (optional)}{
#'     A character-valued data.frame of dimension \eqn{D \times ?} with sequence IDs as rownames.
#'     Expected columns include taxonomic levels such as \code{Kingdom}, \code{Phylum},
#'     \code{Class}, \code{Order}, \code{Genus}, \code{Species}, and optionally \code{Strain}.
#'     Each row describes the taxonomy assigned to a sequence ID. A column named \code{Sequence}
#'     should link sequence IDs to their actual DNA/RNA sequence (e.g., 16S region).
#'
#'     Taxa not classified at a given level should be prefixed with \code{uc_}
#'     followed by the taxonomic level and classification name
#'     (e.g., \code{uc_p_Bacteroidota} for an unclassified phylum).
#'   }
#'
#'   \item{\code{phylo} (optional)}{
#'     A phylogenetic tree object (class to be standardized).
#'     Let the maintainers know if your dataset includes a tree so that
#'     we can determine a common format.
#'   }
#' }
#'
#' @examples
#' ## Run every parser
#' all_data <- microbialscalerepo(base_directory = "data_repository")
#'
#' ## Run two specific parsers and rename outputs
#' subset <- microbialscalerepo(
#'   studies        = c(Vandeputte2021 = "2021_vandeputte_naturecommunications_flow_timeseries",
#'                      Pereira2023   = "2023_pereira_nature_nervous"),
#'   base_directory = "data_repository",
#'   align_samples  = TRUE,
#'   save_to        = "my_parsed.RData"
#' )
#' names(subset)
#' 
#' @export
# ─────────────────────────────────────────────────────────────────────────────
microbialscalerepo <- function(
  studies        = NULL,
  base_directory = ".",
  rawdata        = FALSE,
  align_samples  = FALSE,
  save_to        = NULL
) {
  # --- helper: coerce list → named character vector --------------------------
  if (is.list(studies)) {
    studies <- unlist(studies, use.names = TRUE)
  }

  # ---------------------------- Check
  stopifnot(is.logical(rawdata), length(rawdata) == 1)
  stopifnot(is.logical(align_samples), length(align_samples) == 1)

  # --- locate parser sub‑folders --------------------------------------------
  setwd(base_directory)
  all_dirs    <- list.dirs(base_directory, full.names = FALSE, recursive = FALSE)
  parser_dirs <- all_dirs[file.exists(file.path(base_directory, all_dirs, "parse.R"))]

  # --- resolve selections ----------------------------------------------------
  if (is.null(studies)) {
    selected  <- parser_dirs
    out_names <- parser_dirs
  } else if (is.null(names(studies))) {
    missing <- setdiff(studies, parser_dirs)
    if (length(missing))
      stop("Parser folder(s) not found: ", paste(missing, collapse = ", "))
    selected  <- studies
    out_names <- studies
  } else {
    missing <- setdiff(unname(studies), parser_dirs)
    if (length(missing))
      stop("Parser folder(s) not found: ", paste(missing, collapse = ", "))
    selected  <- unname(studies)
    out_names <- names(studies)
  }

  n <- length(selected)
  if (n == 0) return(invisible(list()))

  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  on.exit(close(pb), add = TRUE)

  parsed_list <- vector("list", n)
  names(parsed_list) <- out_names
  
  # Create a list to store validation results
  validation_results <- vector("list", n)
  names(validation_results) <- out_names

  # --- main loop -------------------------------------------------------------
  for (i in seq_along(selected)) {
    parser          <- selected[i]
    parse_file      <- file.path(base_directory, parser, "parse.R")
    helperfunctions <- file.path(base_directory, "helpers.R")

    env <- new.env()
    if (file.exists(helperfunctions)) {
    sys.source(helperfunctions, envir = env)
    }
    sys.source(parse_file, envir = env)

    fun_name <- paste0("parse_", parser)
    if (!exists(fun_name, envir = env, mode = "function")) {
      warning("Skipping ", parser, ": function ", fun_name, "() not found")
      utils::setTxtProgressBar(pb, i); next
    }
    
    res <- tryCatch(suppressWarnings(get(fun_name, envir = env)(raw = rawdata, align = align_samples)), error = function(e) e)

    if (inherits(res, "error")) {
      warning("Error in parser '", parser, "': ", res$message)
      utils::setTxtProgressBar(pb, i); next
    }

    # Store the original parsed result
    parsed_list[[i]] <- res
    
    # Collect validation result separately
    validation_results[[i]] <- validate_output_structure(res, study_name = parser)
    
    utils::setTxtProgressBar(pb, i)
  }
  utils::setTxtProgressBar(pb, n)

  # --- optional save ---------------------------------------------------------
  if (!is.null(save_to)) {
    dir.create(dirname(save_to), showWarnings = FALSE, recursive = TRUE)
    save(parsed_list, file = save_to)
    
    # Save validation results to a separate file
    validation_file <- sub("\\.RData$", "_validation.RData", save_to)
    save(validation_results, file = validation_file)
    
    # Generate and save validation summary
    summary_file <- sub("\\.RData$", "_validation_summary.txt", save_to)
    sink(summary_file)
    cat("Validation Summary\n")
    cat("=================\n\n")
    
    for (study in names(validation_results)) {
      if (!is.null(validation_results[[study]])) {
        cat(sprintf("Study: %s\n", study))
        cat("Structure:\n")
        for (elem in names(validation_results[[study]])) {
          cat(sprintf("  %s: %s\n", elem, validation_results[[study]][[elem]]))
        }
        cat("\n")
      }
    }
    sink()
  }

  parsed_list
}

# ─────────────────────────────────────────────────────────────────────────────
# Command‑line interface when run via Rscript:
# ─────────────────────────────────────────────────────────────────────────────
if (!interactive()) {
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option(c("-s", "--studies"), type = "character",
                help = "Comma-separated list of parsers (names or paths)"),
    make_option(c("-n", "--names"), type = "character",
                help = "Comma-separated output names (optional, same length)"),
    make_option(c("-b", "--base"),   type = "character", default = ".",
                help = "Base directory [default %default]"),
    make_option(c("-r", "--raw"),    action = "store_true", default = FALSE,
                help = "Pass raw = TRUE to parsers"),
    make_option(c("-a", "--align"),  action = "store_true", default = FALSE,
                help = "Align samples to scale IDs"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Path to save .RData output" )
  )
  parser <- OptionParser(option_list = option_list,
                         description = "Run microbialscalerepo from the command line")
  opts <- parse_args(parser)

  studs <- if (!is.null(opts$studies))
            trimws(strsplit(opts$studies, ",")[[1]])
        else NULL
    if (!is.null(opts$names)) {
        nm <- trimws(strsplit(opts$names, ",")[[1]])
    if (length(nm) != length(studs))
        stop("--names and --studies must have the same number of entries")
    studs <- setNames(studs, nm)
  }

  res <- microbialscalerepo(
    studies        = studs,
    base_directory = opts$base,
    rawdata        = opts$raw,
    align_samples  = opts$align,
    save_to        = opts$output
  )

  if (is.null(opts$output)) {
    cat("Parsed studies:\n", paste(names(res), collapse = "\n"), "\n")
  } else {
    cat("Saved parsed list to:", opts$output, "\n")
}
}