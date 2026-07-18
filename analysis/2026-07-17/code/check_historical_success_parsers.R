#!/usr/bin/env Rscript

library(mutt)

data_root <- Sys.getenv("MUTT_DATA_DIR", unset = "")
output_dir <- Sys.getenv("MUTT_CHECK_OUTPUT", unset = "")
if (!nzchar(data_root) || !dir.exists(data_root)) {
  stop("MUTT_DATA_DIR must name the studies directory.", call. = FALSE)
}
if (!nzchar(output_dir)) stop("MUTT_CHECK_OUTPUT is required.", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

repository_root <- dirname(normalizePath(data_root, mustWork = TRUE))
historical_validation_path <- Sys.getenv(
  "MUTT_HISTORICAL_VALIDATION",
  unset = file.path(repository_root, "datasetsfromrepo_validation.RData")
)
if (!file.exists(historical_validation_path)) {
  stop("Historical validation artifact not found: ", historical_validation_path, call. = FALSE)
}

baseline <- data.frame(
  output_name = c(
    "Vandeputte2021", "CvandeVelde2022", "Vandeputte2017", "Pereira2023",
    "Krawczyk2022", "Liao2021", "Stammler2016", "Dreier2022", "Marotz2021",
    "Vieira_Silva2019", "Contijoch2019", "Tunsakul2024", "Alessandri2024",
    "Maghini2023", "Garcia_Martinez2024", "Sternes2024", "Rao2021",
    "Tettamanti_Boshier2020", "Kruger2024", "Liu2017", "Jin2022",
    "Zaramela2022", "Feng2023", "Reese2022", "Barlow2020", "Morton2019",
    "Prochazkova2024", "Zemb2020", "Galazzo2020", "Lin2019", "Suriano2022",
    "Thiruppathy2025", "Kallastu2023"
  ),
  study = c(
    "2021_vandeputte_naturecommunications_flow_timeseries",
    "2022_cvandevelde_ismecommunications_culturedflowhumanfecal",
    "2017_vandeputte_nature_flow",
    "2023_pereira_nature_nervous",
    "2022_krawczyk_microbiome_tickgeographicaldistributionqpcr",
    "2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct",
    "2016_stammler_microbiome_micehuman",
    "2022_dreier_bmcmicrobiology_cheeseqpcr",
    "2021_marotz_mSystems_oral_mouthwash",
    "2019_vieirasilva_naturemicrobiology_pscibd",
    "2019_contijoch_elife_multispeciesqPCRshotgunandamplicon",
    "2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity",
    "2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota",
    "2023_maghini_naturebiotechnology_samplemesurement",
    "2024_garciamartinez_bmcmicrobiology_ckdanddysbiosiswithserum",
    "2024_sternes_frontmicrobiol_IBDppiqPCR",
    "2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification",
    "2020_tettamantiboshier_msystems_vaginaltimeseries",
    "2024_kruger_scientificreports_ddpcrhealthysubjects",
    "2017_liu_mbio_penilehivqPCR",
    "2022_jin_natureComm_technicalReplicates",
    "2022_zaramela_msystems_synDNA",
    "2023_feng_imetawiley_chickensegment",
    "2021_reese_cell_chimpanzee",
    "2020_barlow_naturecommunications_miceGI",
    "2019_morton_naturecommunications_songbird_oral",
    "2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal",
    "2020_zemb_microOpen_spike",
    "2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy",
    "2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein",
    "2022_suriano_aps_micefecal",
    "2025_thiruppathy_microbiome_relicDNAflow",
    "2023_kallastu_research_foodscience_food"
  ),
  stringsAsFactors = FALSE
)
requested <- trimws(Sys.getenv("MUTT_HISTORICAL_STUDIES", unset = ""))
if (nzchar(requested)) {
  requested <- trimws(strsplit(requested, ",", fixed = TRUE)[[1L]])
  known <- c(baseline$output_name, baseline$study)
  unknown <- setdiff(requested, known)
  if (length(unknown)) {
    stop("Unknown historical studies: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  baseline <- baseline[
    baseline$output_name %in% requested | baseline$study %in% requested,
    ,
    drop = FALSE
  ]
}

historical_env <- new.env(parent = emptyenv())
loaded <- load(historical_validation_path, envir = historical_env)
if (!"validation_results" %in% loaded) {
  stop("Historical validation artifact does not contain validation_results.", call. = FALSE)
}
historical <- historical_env$validation_results
missing_baselines <- setdiff(baseline$output_name, names(historical))
if (length(missing_baselines)) {
  stop("Missing historical outputs: ", paste(missing_baselines, collapse = ", "), call. = FALSE)
}

core_elements <- c("counts", "proportions", "tax", "scale", "metadata")
normalize_tax_for_historical_validation <- function(x, expected_structure) {
  tokens <- regmatches(
    expected_structure,
    gregexpr("DATAFRAME\\([0-9]+ x [0-9]+\\)", expected_structure)
  )[[1L]]
  expected_rows <- as.integer(sub("DATAFRAME\\(([0-9]+) x .*", "\\1", tokens))
  expected_cols <- as.integer(sub(".* x ([0-9]+)\\)", "\\1", tokens))
  position <- 0L
  walk <- function(value) {
    if (is.data.frame(value)) {
      position <<- position + 1L
      if (
        position <= length(expected_cols) &&
        nrow(value) == expected_rows[[position]] &&
        ncol(value) == expected_cols[[position]] + 1L &&
        "Sequence" %in% names(value)
      ) {
        value <- value[, setdiff(names(value), "Sequence"), drop = FALSE]
      }
      return(value)
    }
    if (is.list(value)) return(lapply(value, walk))
    value
  }
  walk(x)
}
baseline$core_structure <- vapply(
  baseline$output_name,
  function(name) paste(unlist(historical[[name]][core_elements]), collapse = " || "),
  character(1)
)
baseline$align_samples <- TRUE
baseline$source_artifact <- basename(historical_validation_path)
utils::write.csv(
  baseline,
  file.path(output_dir, "historical_success_baseline.csv"),
  row.names = FALSE,
  na = ""
)

report_rows <- vector("list", nrow(baseline))
for (i in seq_len(nrow(baseline))) {
  output_name <- baseline$output_name[[i]]
  study <- baseline$study[[i]]
  cat(sprintf("[%d/%d] %s\n", i, nrow(baseline), study))
  emitted <- character()
  started <- Sys.time()
  result <- withCallingHandlers(
    tryCatch(
      mutt(
        studies = stats::setNames(study, output_name),
        base_directory = data_root,
        align_samples = TRUE,
        functional = FALSE
      ),
      error = identity
    ),
    warning = function(w) {
      emitted <<- c(emitted, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  status <- "error"
  structure_match <- FALSE
  mismatch <- ""
  parser_error <- ""
  parser_warnings <- ""
  data_source <- ""
  runtime_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  if (inherits(result, "error")) {
    parser_error <- conditionMessage(result)
  } else {
    audit <- attr(result, "audit")
    status <- audit$status[[1L]]
    data_source <- audit$data_source[[1L]]
    runtime_seconds <- audit$runtime_seconds[[1L]]
    parser_error <- audit$error[[1L]]
    parser_warnings <- audit$warnings[[1L]]
    current <- attr(result, "validation")[[output_name]][core_elements]
    parser_tax <- normalize_tax_for_historical_validation(
      result[[output_name]]$tax,
      historical[[output_name]]$tax
    )
    current$tax <- mutt:::validate_output_structure(list(tax = parser_tax))$tax
    expected <- historical[[output_name]][core_elements]
    comparison <- all.equal(current, expected, check.attributes = TRUE)
    structure_match <- isTRUE(comparison)
    if (!structure_match) mismatch <- paste(comparison, collapse = " | ")
  }

  warning_messages <- unique(c(emitted, parser_warnings))
  warning_messages <- warning_messages[!is.na(warning_messages) & nzchar(warning_messages)]
  report_rows[[i]] <- data.frame(
    output_name = output_name,
    study = study,
    status = status,
    align_samples = TRUE,
    core_structure_match = structure_match,
    runtime_seconds = runtime_seconds,
    data_source = data_source,
    warnings = paste(warning_messages, collapse = " | "),
    error = parser_error,
    mismatch = mismatch,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %s; structure_match=%s (%.1f s)\n", status, structure_match, runtime_seconds))
  rm(result)
  invisible(gc())
}

report <- do.call(rbind, report_rows)
utils::write.csv(
  report,
  file.path(output_dir, "historical_success_check.csv"),
  row.names = FALSE,
  na = ""
)
summary <- data.frame(
  historical_studies = nrow(report),
  parser_successes = sum(report$status == "success"),
  structure_matches = sum(report$core_structure_match),
  failures = sum(report$status != "success" | !report$core_structure_match)
)
utils::write.csv(
  summary,
  file.path(output_dir, "historical_success_check_summary.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "historical_success_session_info.txt")
)
print(summary, row.names = FALSE)
if (summary$failures[[1L]] > 0L) quit(save = "no", status = 1L)
