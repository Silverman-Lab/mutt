# ─────────────────────────────────────────────────────────────────────────────
#' Microbial Scale Repository
#'
#' Public genomic (Amplicon or Shotgun Sequencing) datasets parsed or reanalyzed
#' stored in subdirectories under `base_directory` that contain a `parse.R`
#' file, source each in isolation, run its `parse_<folder>()` function, and
#' collect the results. 
#'
#' @author Maxwell Konnaris
#' @author Justin Silverman
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
#' @param verbose Logical(1). If `TRUE`, verbose output to console of the structure
#'        of returned studies named list
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
#' 
#'   \item{\code{studydemographics} (optional)}{
#'     List. Manually curated information about a particular studies data.
#'   }
#' }
#'
#' @examples
#' ## Run every parser
#' all_data <- totallia(base_directory = "totallia")
#'
#' ## Run two specific parsers and rename outputs
#' subset <- mutt(
#'   studies        = c(Vandeputte2021 = "2021_vandeputte_naturecommunications_flow_timeseries",
#'                      Pereira2023   = "2023_pereira_nature_nervous"),
#'   base_directory = "totallia",
#'   align_samples  = TRUE,
#'   save_to        = "my_parsed.RData",
#'   verbose        = TRUE
#' )
#' names(subset)
#' 
#' @export
# ─────────────────────────────────────────────────────────────────────────────
mutt <- function(
  studies        = NULL,
  base_directory = ".",
  rawdata        = FALSE,
  align_samples  = FALSE,
  save_to        = NULL,
  verbose        = FALSE
) {
  # --- helper: coerce list → named character vector --------------------------
  if (is.list(studies)) {
    studies <- unlist(studies, use.names = TRUE)
  }

  ## ---------- helper: recursively drop all‑NA columns -----------------------
  remove_all_na <- function(x) {
    if (is.data.frame(x)) {
      keep <- vapply(x, function(col) !all(is.na(col)), logical(1))
      x[, keep, drop = FALSE]
    } else if (is.list(x)) {
      lapply(x, remove_all_na)
    } else {
      x
    }
  }

  # ---------------------------- Check ----------------------------------------
  stopifnot(is.logical(rawdata), length(rawdata) == 1)
  stopifnot(is.logical(align_samples), length(align_samples) == 1)

  # --- locate parser sub‑folders ---------------------------------------------
  setwd(base_directory)
  all_dirs    <- list.dirs(base_directory, full.names = FALSE, recursive = FALSE)
  parser_dirs <- all_dirs[file.exists(file.path(base_directory, all_dirs, "parse.R"))]
  if (!file.exists("R/demographics.R")) {
    stop("R/demographics.R not found in base_directory: ", base_directory)
  }
  source("R/demographics.R") 

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
  microbialscalerepository <- vector("list", n)
  names(microbialscalerepository) <- out_names
  validation_results <- vector("list", n)
  names(validation_results) <- out_names

  # --- main loop -------------------------------------------------------------
  for (i in seq_along(selected)) {
    parser          <- selected[i]
    parse_file      <- file.path(base_directory, parser, "parse.R")
    helperfunctions <- file.path(base_directory, "R/helpers.R")

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

    if (verbose) {
      cat(sprintf("\r[%d/%d] %s", i, n, parser))
      flush.console()  
    }

    res <- tryCatch(
      suppressWarnings(suppressMessages(get(fun_name, envir = env)(raw = rawdata, align = align_samples))), 
      error = function(e) {
        warning(sprintf("\n\n----------------------------------------\nError in parser '%s': %s\n----------------------------------------\n", 
                       parser, 
                       e$message))
        e
      }
    )

    if(verbose) {
      if (inherits(res, "error")) {
        warning(sprintf("\n\n----------------------------------------\nError in parser '%s': %s\n----------------------------------------\n", 
                       parser, 
                       res$message))
        utils::setTxtProgressBar(pb, i); next
      }
    }

    res <- remove_all_na(res)
    res <- standardize_output_order(res)
    if (!is.null(datasets[[parser]])) {
        res$studydemographics <- datasets[[parser]]
      }
    microbialscalerepository[[i]] <- res
    validation_results[[i]] <- validate_output_structure(res, study_name = parser)
    
    utils::setTxtProgressBar(pb, i)
  }
  utils::setTxtProgressBar(pb, n)
  close(pb)
  if (verbose) cat("\n")

  microbialscalerepository <-lapply(microbialscalerepository, function(study) {
    if (!is.null(study$tax)) {
      study$tax <- add_sequence_column(study$tax)
    }
    study
  })

  # --- optional save ---------------------------------------------------------
  if (!is.null(save_to)) {
    dir.create(dirname(save_to), showWarnings = FALSE, recursive = TRUE)
    save(microbialscalerepository, file = save_to)
    cat(sprintf("File saved to: %s\n", save_to))
    validation_file <- sub("\\.RData$", "_validation.RData", save_to)
    save(validation_results, file = validation_file)
    if (verbose) {
      summary_file <- sub("\\.RData$", "_validation_summary.txt", save_to)
      summary_text <- capture.output({
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
      })
      writeLines(summary_text, summary_file)
      cat(summary_text, sep="\n")
    }
  }
  microbialscalerepository
}
