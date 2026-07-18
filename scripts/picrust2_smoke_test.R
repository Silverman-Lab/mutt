#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
package_root <- dirname(dirname(script_path))

parser <- optparse::OptionParser(
  description = "Run a small, real, verbose PICRUSt2 smoke test from a MUTT study.",
  option_list = list(
    optparse::make_option(
      "--study-dir", type = "character", dest = "study_dir",
      help = "Study directory containing parse.R and ASV input files."
    ),
    optparse::make_option(
      "--source-id", type = "character", dest = "source_id", default = NULL,
      help = "Optional ASV source ID. The first eligible source is used by default."
    ),
    optparse::make_option(
      "--samples", type = "integer", default = 3L,
      help = "Number of highest-depth samples to use [default: %default]."
    ),
    optparse::make_option(
      "--asvs", type = "integer", default = 3L,
      help = "Number of most abundant ASVs in the selected samples [default: %default]."
    ),
    optparse::make_option(
      "--processes", type = "integer", default = 8L,
      help = "PICRUSt2 worker processes [default: %default]."
    ),
    optparse::make_option(
      "--output-dir", type = "character", dest = "output_dir", default = NULL,
      help = "New output directory. Defaults to study/functional/smoke_test/<timestamp>."
    ),
    optparse::make_option(
      "--prepare-only", action = "store_true", dest = "prepare_only", default = FALSE,
      help = "Prepare and validate inputs, print the command, and stop before PICRUSt2."
    )
  )
)
opt <- optparse::parse_args(parser)

if (is.null(opt$study_dir) || !nzchar(opt$study_dir)) {
  optparse::print_help(parser)
  stop("`--study-dir` is required.", call. = FALSE)
}
for (name in c("samples", "asvs", "processes")) {
  value <- opt[[name]]
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop("`--", name, "` must be a positive integer.", call. = FALSE)
  }
}

study_dir <- normalizePath(opt$study_dir, mustWork = TRUE)
if (!file.exists(file.path(study_dir, "parse.R"))) {
  stop("The study directory does not contain parse.R: ", study_dir, call. = FALSE)
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("This source-tree smoke test requires the pkgload R package.", call. = FALSE)
}
pkgload::load_all(package_root, quiet = TRUE)
ns <- asNamespace("mutt")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- if (is.null(opt$output_dir)) {
  file.path(study_dir, "functional", "smoke_test", timestamp)
} else {
  path.expand(opt$output_dir)
}
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Refusing to overwrite existing output: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

timings <- data.frame(phase = character(), status = character(), seconds = numeric())
write_timings <- function() {
  utils::write.table(
    timings, file.path(output_dir, "timings.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
run_phase <- function(label, code) {
  cat("\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] START: ", label, "\n", sep = "")
  started <- proc.time()[["elapsed"]]
  tryCatch(
    {
      value <- force(code)
      elapsed <- proc.time()[["elapsed"]] - started
      timings <<- rbind(
        timings,
        data.frame(phase = label, status = "passed", seconds = elapsed)
      )
      write_timings()
      cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] DONE:  ", label,
          " (", sprintf("%.1f", elapsed), " s)\n", sep = "")
      value
    },
    error = function(e) {
      elapsed <- proc.time()[["elapsed"]] - started
      timings <<- rbind(
        timings,
        data.frame(phase = label, status = "failed", seconds = elapsed)
      )
      write_timings()
      cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] FAILED: ", label,
          " (", sprintf("%.1f", elapsed), " s)\n", sep = "")
      stop(e)
    }
  )
}

cat("MUTT PICRUSt2 real smoke test\n")
cat("Package root: ", package_root, "\n", sep = "")
cat("Study:       ", study_dir, "\n", sep = "")
cat("Output:      ", output_dir, "\n", sep = "")
cat("Requested:   ", opt$samples, " samples, ", opt$asvs, " ASVs, ",
    opt$processes, " processes\n", sep = "")

prepared <- run_phase("discover, validate, and subset real ASV inputs", {
  sources <- get(".discover_picrust_sources", ns)(study_dir)
  if (!length(sources)) stop("No ASV sequence source was found.", call. = FALSE)
  ids <- vapply(sources, `[[`, character(1), "id")
  source_index <- if (is.null(opt$source_id)) 1L else match(opt$source_id, ids)
  if (is.na(source_index)) {
    stop("Unknown `--source-id`. Available sources: ", paste(ids, collapse = ", "), call. = FALSE)
  }
  source <- sources[[source_index]]
  loaded <- get(".load_picrust_source", ns)(source)
  if (!isTRUE(loaded$prokaryotic)) stop("Selected source is not prokaryotic 16S data.", call. = FALSE)

  sample_order <- names(sort(rowSums(loaded$counts), decreasing = TRUE))
  sample_ids <- head(sample_order, min(opt$samples, length(sample_order)))
  candidate_totals <- colSums(loaded$counts[sample_ids, , drop = FALSE])
  candidate_totals <- sort(candidate_totals[candidate_totals > 0], decreasing = TRUE)
  domain_column <- names(loaded$taxonomy)[
    tolower(names(loaded$taxonomy)) %in% c("kingdom", "domain", "superkingdom")
  ][1L]
  if (!length(domain_column) || is.na(domain_column)) {
    stop("Taxonomy has no Kingdom, Domain, or Superkingdom column.", call. = FALSE)
  }
  domains <- tolower(trimws(as.character(loaded$taxonomy[[domain_column]])))
  classified_ids <- loaded$taxonomy$original_feature_id[
    domains %in% c("bacteria", "archaea", "k__bacteria", "k__archaea",
                   "d__bacteria", "d__archaea")
  ]
  candidate_totals <- candidate_totals[names(candidate_totals) %in% classified_ids]
  asv_ids <- head(names(candidate_totals), min(opt$asvs, length(candidate_totals)))
  if (!length(asv_ids)) {
    stop("Selected samples contain no nonzero classified bacterial/archaeal ASVs.", call. = FALSE)
  }

  counts <- loaded$counts[sample_ids, asv_ids, drop = FALSE]
  sequences <- loaded$sequences[asv_ids]
  taxonomy <- loaded$taxonomy[
    match(asv_ids, loaded$taxonomy$original_feature_id), , drop = FALSE
  ]
  if (anyNA(taxonomy$original_feature_id) ||
      !identical(as.character(taxonomy$original_feature_id), asv_ids)) {
    stop("Smoke-test ASV IDs do not align with taxonomy IDs.", call. = FALSE)
  }
  list(
    source = source,
    counts = counts,
    sequences = sequences,
    taxonomy = taxonomy
  )
})

tools <- run_phase("discover and verify PICRUSt2 dependencies", {
  tools <- get(".discover_functional_tools", ns)()
  required <- c("place_seqs", "hsp", "metagenome", "pathway", "biom", "hmmalign")
  missing <- required[
    !nzchar(tools$picrust2_dependencies[required]) |
      is.na(tools$picrust2_dependencies[required])
  ]
  if (!nzchar(tools$picrust2) || length(missing)) {
    stop(
      "PICRUSt2 is unavailable or incomplete. Missing: ",
      paste(c(if (!nzchar(tools$picrust2)) "picrust2_pipeline.py", missing), collapse = ", "),
      call. = FALSE
    )
  }
  cat("PICRUSt2: ", tools$picrust2_version, "\n", sep = "")
  cat("BIOM:     ", tools$biom_version, "\n", sep = "")
  tools
})

input_dir <- file.path(output_dir, "input")
orientation <- run_phase("compare supplied and reverse-complement ASV orientations", {
  get(".orient_picrust_sequences", ns)(
    prepared$sequences,
    file.path(input_dir, "orientation"),
    tools
  )
})
if (!any(orientation$audit$passes_min_align)) {
  stop("No selected ASV passed the 0.8 alignment threshold.", call. = FALSE)
}

inputs <- run_phase("write FASTA and validated HDF5 BIOM inputs", {
  inputs <- get(".write_picrust_inputs", ns)(
    prepared$counts,
    orientation$sequences,
    input_dir,
    tools$picrust2_dependencies[["biom"]]
  )
  rows <- match(inputs$mapping_data$original_feature_id, orientation$audit$original_feature_id)
  inputs$mapping_data <- cbind(
    inputs$mapping_data,
    orientation$audit[rows, setdiff(names(orientation$audit), "original_feature_id"), drop = FALSE]
  )
  inputs$mapping_data$taxonomy_row_id <- inputs$mapping_data$original_feature_id
  utils::write.table(
    inputs$mapping_data, inputs$mapping, sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    prepared$taxonomy, file.path(input_dir, "asv_taxonomy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  inputs
})

raw_dir <- file.path(output_dir, "raw")
picrust_args <- c(
  "-s", inputs$fasta,
  "-i", inputs$table,
  "-o", raw_dir,
  "-p", as.character(opt$processes),
  "--min_align", "0.8",
  "--coverage",
  "--stratified",
  "--verbose"
)
command_text <- paste(c(shQuote(tools$picrust2), vapply(picrust_args, shQuote, character(1))), collapse = " ")
writeLines(command_text, file.path(output_dir, "command.txt"))
writeLines(
  c(
    R.version.string,
    paste("mutt", as.character(utils::packageVersion("mutt"))),
    tools$picrust2_version,
    tools$biom_version
  ),
  file.path(output_dir, "session_info.txt")
)

cat("\nExact PICRUSt2 command:\n", command_text, "\n", sep = "")
if (isTRUE(opt$prepare_only)) {
  cat("\nPREPARATION PASSED; PICRUSt2 was not run because --prepare-only was set.\n")
  cat("Output directory: ", output_dir, "\n", sep = "")
  quit(save = "no", status = 0L)
}
cat("\nPICRUSt2 will now print its wrapped placement, HSP, metagenome, and pathway commands.\n")
cat("The EPA-NG reference placement has a large fixed memory/runtime cost even for three ASVs.\n")
Sys.setenv(PYTHONUNBUFFERED = "1")
run_args <- picrust_args
run_args[c(2L, 4L, 6L)] <- vapply(run_args[c(2L, 4L, 6L)], shQuote, character(1))
run_phase("run the full real PICRUSt2 pipeline", {
  status <- system2(tools$picrust2, run_args, stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) {
    stop("PICRUSt2 exited with status ", status, ". Partial files remain in ", output_dir,
         call. = FALSE)
  }
  invisible(status)
})

result <- run_phase("parse and validate all required outputs", {
  result <- get(".parse_picrust_output", ns)(raw_dir, prepared$counts)
  result$asv_mapping <- inputs$mapping_data
  result$taxonomy <- prepared$taxonomy
  result$provenance <- list(
    source = prepared$source$counts_file,
    tool_version = tools$picrust2_version,
    output_directory = output_dir,
    min_align = 0.8,
    processes = opt$processes,
    smoke_test = TRUE
  )

  read_contributions <- get("read_picrust2_contributions", ns)
  for (type in c("ec", "ko", "metacyc_abundance")) {
    preview <- read_contributions(result, type, n_max = 10L)
    if (!nrow(preview) || anyNA(preview$original_feature_id) ||
        !all(preview$original_feature_id %in% colnames(prepared$counts))) {
      stop("Invalid stratified contribution mapping for ", type, ".", call. = FALSE)
    }
  }
  result
})

saveRDS(result, file.path(output_dir, "result.rds"))
write_timings()

cat("\nSMOKE TEST PASSED\n")
cat("Samples:            ", nrow(prepared$counts), "\n", sep = "")
cat("ASVs:               ", ncol(prepared$counts), "\n", sep = "")
cat("EC features:        ", result$qc$ec_features, "\n", sep = "")
cat("KO features:        ", result$qc$ko_features, "\n", sep = "")
cat("MetaCyc pathways:   ", result$qc$metacyc_pathways, "\n", sep = "")
cat("Output directory:   ", output_dir, "\n", sep = "")
cat("Command record:     ", file.path(output_dir, "command.txt"), "\n", sep = "")
cat("Timing record:      ", file.path(output_dir, "timings.tsv"), "\n", sep = "")
