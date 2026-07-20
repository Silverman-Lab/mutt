.functional_mode <- function(functional) {
  if (isFALSE(functional)) return("off")
  if (isTRUE(functional)) return("use")
  if (is.character(functional) && length(functional) == 1L &&
      !is.na(functional) && tolower(functional) == "rebuild") {
    return("rebuild")
  }
  stop("`functional` must be FALSE, TRUE, or \"rebuild\".", call. = FALSE)
}

.functional_engine_version <- "5"
.picrust_min_align <- 0.8

.functional_branch_modality <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return("unknown")
  tokens <- unlist(
    strsplit(tolower(gsub("\\\\", "/", path)), "[^a-z0-9]+", perl = TRUE),
    use.names = FALSE
  )
  tokens <- tokens[nzchar(tokens)]
  metagenomic <- c(
    "shotgun", "metagenome", "metagenomic", "metaphlan", "metaphlan2",
    "metaphlan3", "metaphlan4", "motu", "motus", "kraken", "kraken2",
    "bracken", "humann", "humann2", "humann3", "wgs"
  )
  amplicon <- c(
    "amplicon", "dada2", "asv", "asvs", "otu", "otus", "rdp", "rdp16",
    "rdp19", "16s", "18s", "its"
  )
  if (any(tokens %in% metagenomic)) return("metagenomic")
  if (any(tokens %in% amplicon)) return("amplicon")
  "unknown"
}

.functional_faprotax_exclusion_template <- function() {
  data.frame(
    branch = character(),
    reason = character(),
    source = character(),
    samples = integer(),
    features = integer(),
    stringsAsFactors = FALSE
  )
}

.functional_processes <- function() {
  requested <- suppressWarnings(as.integer(Sys.getenv("MUTT_FUNCTIONAL_PROCESSES", "")))
  if (length(requested) == 1L && !is.na(requested) && requested > 0L) {
    return(requested)
  }
  detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (length(detected) != 1L || is.na(detected) || detected < 1L) detected <- 1L
  max(1L, min(8L, as.integer(detected)))
}

.functional_manifest_template <- function() {
  data.frame(
    method = character(),
    branch = character(),
    status = character(),
    reason = character(),
    source = character(),
    samples = integer(),
    features = integer(),
    input_hash = character(),
    tool_version = character(),
    runtime_seconds = numeric(),
    output_directory = character(),
    stringsAsFactors = FALSE
  )
}

.empty_functional_result <- function(requested = FALSE) {
  out <- list(
    picrust2 = list(),
    faprotax = list(),
    manifest = .functional_manifest_template()
  )
  attr(out, "requested") <- isTRUE(requested)
  out
}

.manifest_row <- function(method, branch, status, reason = "", source = "",
                          samples = NA_integer_, features = NA_integer_,
                          input_hash = "", tool_version = "",
                          runtime_seconds = NA_real_, output_directory = "") {
  data.frame(
    method = method,
    branch = branch,
    status = status,
    reason = reason,
    source = source,
    samples = as.integer(samples),
    features = as.integer(features),
    input_hash = input_hash,
    tool_version = tool_version,
    runtime_seconds = as.numeric(runtime_seconds),
    output_directory = output_directory,
    stringsAsFactors = FALSE
  )
}

.safe_branch_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", sub("_+$", "", x))
  if (!nzchar(x)) "dataset" else x
}

.set_nested <- function(x, path, value) {
  path <- path[nzchar(path)]
  if (!length(path)) return(value)
  key <- path[[1L]]
  if (length(path) == 1L) {
    x[[key]] <- value
  } else {
    child <- x[[key]]
    if (!is.list(child)) child <- list()
    x[[key]] <- .set_nested(child, path[-1L], value)
  }
  x
}

.flatten_tables <- function(x, path = character()) {
  if (is.data.frame(x) || is.matrix(x)) {
    if (!length(dim(x)) || any(dim(x) == 0L)) return(list())
    return(setNames(list(x), paste(path, collapse = "/")))
  }
  if (!is.list(x) || !length(x)) return(list())
  out <- list()
  nms <- names(x)
  if (is.null(nms)) nms <- as.character(seq_along(x))
  for (i in seq_along(x)) {
    out <- c(out, .flatten_tables(x[[i]], c(path, nms[[i]])))
  }
  out
}

.object_md5 <- function(x) {
  path <- tempfile("mutt_hash_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

.files_md5 <- function(paths) {
  paths <- normalizePath(paths[file.exists(paths)], mustWork = TRUE)
  if (!length(paths)) stop("No source files were available for hashing.", call. = FALSE)
  .object_md5(as.list(unname(tools::md5sum(paths))))
}

.command_invocation <- function(command, args = character()) {
  if (nzchar(command) && grepl("\\.py$", command, ignore.case = TRUE) &&
      file.exists(command) && file.access(command, mode = 1L) != 0L) {
    python <- unname(Sys.which("python3"))
    if (!nzchar(python)) python <- unname(Sys.which("python"))
    if (!nzchar(python)) {
      stop("Python was not found on PATH for script: ", command, call. = FALSE)
    }
    return(list(command = python, args = c(shQuote(command), args)))
  }
  list(command = command, args = args)
}

.command_version <- function(command) {
  if (!nzchar(command)) return("")
  for (flag in c("--version", "-v")) {
    invocation <- tryCatch(.command_invocation(command, flag), error = function(e) NULL)
    if (is.null(invocation)) next
    out <- tryCatch(
      suppressWarnings(system2(invocation$command, invocation$args, stdout = TRUE, stderr = TRUE)),
      error = function(e) character()
    )
    status <- attr(out, "status", exact = TRUE)
    if (length(out) && (is.null(status) || identical(as.integer(status), 0L))) {
      return(trimws(out[[1L]]))
    }
  }
  "unknown"
}

.discover_picrust_reference_files <- function(picrust2) {
  empty <- setNames(rep("", 4L), c("bacteria_msa", "bacteria_hmm", "archaea_msa", "archaea_hmm"))
  if (!nzchar(picrust2)) return(empty)
  python <- unname(Sys.which("python"))
  if (!nzchar(python)) python <- unname(Sys.which("python3"))
  if (!nzchar(python)) return(empty)

  code <- paste(
    "import os, picrust2",
    "root = os.path.join(os.path.dirname(picrust2.__file__), 'default_files')",
    "print(os.path.join(root, 'bacteria', 'bac_ref', 'bac_ref.fna'))",
    "print(os.path.join(root, 'bacteria', 'bac_ref', 'bac_ref.hmm'))",
    "print(os.path.join(root, 'archaea', 'arc_ref', 'arc_ref.fna'))",
    "print(os.path.join(root, 'archaea', 'arc_ref', 'arc_ref.hmm'))",
    sep = "; "
  )
  paths <- tryCatch(
    suppressWarnings(system2(python, c("-c", shQuote(code)), stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (length(paths) != 4L || any(!file.exists(paths))) return(empty)
  setNames(normalizePath(paths), names(empty))
}

.discover_functional_tools <- function() {
  picrust2 <- unname(Sys.which("picrust2_pipeline.py"))
  picrust_commands <- c(
    place_seqs = "place_seqs.py",
    hsp = "hsp.py",
    metagenome = "metagenome_pipeline.py",
    pathway = "pathway_pipeline.py",
    biom = "biom",
    hmmalign = "hmmalign"
  )
  picrust_dependencies <- unname(Sys.which(picrust_commands))
  names(picrust_dependencies) <- names(picrust_commands)

  faprotax_dir <- Sys.getenv("FAPROTAX_DIR", unset = "")
  bundled_faprotax_dirs <- c(
    system.file("extdata", "FAPROTAX_1.2.12", package = "mutt"),
    system.file("helperdata", "FAPROTAX_1.2.12", package = "mutt"),
    system.file("FAPROTAX_1.2.12", package = "mutt"),
    file.path(getwd(), "helperdata", "FAPROTAX_1.2.12"),
    file.path(getwd(), "inst", "extdata", "FAPROTAX_1.2.12")
  )
  bundled_faprotax_dirs <- unique(
    bundled_faprotax_dirs[nzchar(bundled_faprotax_dirs) & dir.exists(bundled_faprotax_dirs)]
  )
  faprotax <- unname(Sys.which("collapse_table.py"))
  if (!nzchar(faprotax)) {
    script_candidates <- c(
      if (nzchar(faprotax_dir)) file.path(faprotax_dir, "collapse_table.py"),
      file.path(bundled_faprotax_dirs, "collapse_table.py")
    )
    script_candidates <- script_candidates[file.exists(script_candidates)]
    if (length(script_candidates)) faprotax <- normalizePath(script_candidates[[1L]])
  }

  prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  db_candidates <- c(
    if (nzchar(faprotax_dir)) file.path(faprotax_dir, "FAPROTAX.txt"),
    file.path(bundled_faprotax_dirs, "FAPROTAX.txt"),
    if (nzchar(faprotax)) file.path(dirname(faprotax), "FAPROTAX.txt"),
    if (nzchar(faprotax)) file.path(dirname(dirname(faprotax)), "share", "faprotax", "FAPROTAX.txt"),
    if (nzchar(prefix)) file.path(prefix, "share", "faprotax", "FAPROTAX.txt"),
    if (nzchar(prefix)) file.path(prefix, "opt", "faprotax", "FAPROTAX.txt")
  )
  db_candidates <- unique(db_candidates[nzchar(db_candidates)])
  faprotax_db <- db_candidates[file.exists(db_candidates)][1L]
  if (!length(faprotax_db) || is.na(faprotax_db)) faprotax_db <- ""
  faprotax_db_version <- if (nzchar(faprotax_db)) {
    unname(tools::md5sum(faprotax_db))
  } else {
    ""
  }

  list(
    picrust2 = picrust2,
    picrust2_version = .command_version(picrust2),
    picrust2_dependencies = picrust_dependencies,
    picrust2_references = .discover_picrust_reference_files(picrust2),
    biom_version = .command_version(picrust_dependencies[["biom"]]),
    faprotax = faprotax,
    faprotax_db = faprotax_db,
    faprotax_db_version = faprotax_db_version,
    faprotax_version = if (nzchar(faprotax)) .command_version(faprotax) else ""
  )
}

.read_rds_zip <- function(path, pattern) {
  td <- tempfile("mutt_unzip_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  extracted <- utils::unzip(path, exdir = td)
  selected <- extracted[grepl(pattern, basename(extracted), ignore.case = TRUE)]
  if (!length(selected)) stop("Expected file not found in archive: ", path, call. = FALSE)
  readRDS(selected[[1L]])
}

.dna_fraction <- function(x) {
  if (is.null(x) || !length(x)) return(0)
  mean(grepl("^[ACGTRYSWKMBDHVN]+$", x, ignore.case = TRUE))
}

.reverse_complement_dna <- function(x) {
  x <- toupper(as.character(x))
  vapply(
    strsplit(x, "", fixed = TRUE),
    function(z) {
      paste(chartr("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN", rev(z)), collapse = "")
    },
    character(1)
  )
}

.assert_no_reverse_complement_duplicates <- function(sequences, context) {
  sequences <- toupper(as.character(sequences))
  reverse_match <- match(.reverse_complement_dna(sequences), sequences, nomatch = 0L)
  reverse_pairs <- which(reverse_match > seq_along(reverse_match))
  if (length(reverse_pairs) > 0) {
    stop(
      context, " contains ", length(reverse_pairs),
      " exact reverse-complement feature pair(s). Normalize amplicon orientation ",
      "and rebuild counts and taxonomy before PICRUSt2 inference.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.collapse_reverse_complement_duplicates <- function(counts, sequences, tax = NULL) {
  sequences <- toupper(as.character(sequences))
  names(sequences) <- colnames(counts)
  reverse <- .reverse_complement_dna(sequences)
  reverse_match <- match(reverse, sequences, nomatch = 0L)
  paired <- reverse_match != 0L & reverse_match != seq_along(sequences)
  if (!any(paired)) {
    return(list(counts = counts, sequences = sequences, taxonomy = tax, merged_pairs = 0L))
  }

  feature_group <- sequences
  feature_group[paired] <- pmin(sequences[paired], reverse[paired])
  group_order <- unique(feature_group)
  before <- rowSums(counts)
  collapsed <- t(rowsum(t(counts), group = feature_group, reorder = FALSE))
  collapsed <- collapsed[, group_order, drop = FALSE]
  if (!identical(as.numeric(rowSums(collapsed)), as.numeric(before))) {
    stop("Reverse-complement merging changed sample read totals.", call. = FALSE)
  }

  aligned_taxonomy <- .match_taxonomy_to_counts(counts, tax)
  collapsed_taxonomy <- if (is.null(aligned_taxonomy)) {
    NULL
  } else {
    .consensus_taxonomy(aligned_taxonomy, feature_group, group_order)
  }
  collapsed_sequences <- setNames(group_order, group_order)
  .assert_no_reverse_complement_duplicates(
    collapsed_sequences,
    "Reverse-complement-normalized PICRUSt2 sequence input"
  )
  list(
    counts = collapsed,
    sequences = collapsed_sequences,
    taxonomy = collapsed_taxonomy,
    merged_pairs = as.integer(sum(table(feature_group) == 2L))
  )
}

.normalize_count_matrix <- function(x, require_sequences = FALSE) {
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("Count input must be a matrix or data.frame.", call. = FALSE)
  }
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("Count input contains missing, non-finite, or negative values.", call. = FALSE)
  }
  if (any(abs(x - round(x)) > 1e-8)) {
    stop("Functional inference requires integer-valued count data.", call. = FALSE)
  }

  if (require_sequences && .dna_fraction(rownames(x)) > .dna_fraction(colnames(x))) {
    x <- t(x)
  }
  if (require_sequences && .dna_fraction(colnames(x)) < 1) {
    stop("PICRUSt2 input features must have matching DNA sequences.", call. = FALSE)
  }
  if (is.null(rownames(x)) || any(!nzchar(rownames(x))) || anyDuplicated(rownames(x))) {
    stop("Count input requires unique, non-empty sample IDs.", call. = FALSE)
  }
  if (is.null(colnames(x)) || any(!nzchar(colnames(x))) || anyDuplicated(colnames(x))) {
    stop("Count input requires unique, non-empty feature IDs.", call. = FALSE)
  }

  keep <- colSums(x) > 0
  x <- x[, keep, drop = FALSE]
  if (!nrow(x) || !ncol(x)) stop("Count input is empty after removing zero-count features.", call. = FALSE)
  x
}

.normalize_abundance_matrix <- function(x) {
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("Abundance input must be a matrix or data.frame.", call. = FALSE)
  }
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "numeric")
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("Abundance input contains missing, non-finite, or negative values.", call. = FALSE)
  }
  if (is.null(rownames(x)) || any(!nzchar(rownames(x))) || anyDuplicated(rownames(x))) {
    stop("Abundance input requires unique, non-empty sample IDs.", call. = FALSE)
  }
  if (is.null(colnames(x)) || any(!nzchar(colnames(x))) || anyDuplicated(colnames(x))) {
    stop("Abundance input requires unique, non-empty feature IDs.", call. = FALSE)
  }
  keep <- colSums(x) > 0
  x <- x[, keep, drop = FALSE]
  if (!nrow(x) || !ncol(x)) {
    stop("Abundance input is empty after removing zero-abundance features.", call. = FALSE)
  }
  x
}

.rank_columns <- function(tax) {
  wanted <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  nms <- names(tax)
  idx <- match(tolower(wanted), tolower(nms))
  setNames(nms[idx[!is.na(idx)]], wanted[!is.na(idx)])
}

.clean_tax_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "unknown", "Unknown", "unclassified", "Unclassified")] <- NA_character_
  x[grepl("^uc_[kpcofgs]_", x, ignore.case = TRUE)] <- NA_character_
  x <- sub("^[dkpcofgst]__?", "", x, ignore.case = TRUE)
  x <- trimws(x)
  x[!nzchar(x)] <- NA_character_
  x
}

.is_prokaryotic_taxonomy <- function(tax) {
  if (!is.data.frame(tax) && !is.matrix(tax)) return(FALSE)
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)
  ranks <- .rank_columns(tax)
  domain_col <- ranks[c("Kingdom")]
  if (!length(domain_col) || is.na(domain_col)) {
    candidates <- names(tax)[tolower(names(tax)) %in% c("domain", "superkingdom")]
    if (!length(candidates)) return(FALSE)
    domain_col <- candidates[[1L]]
  }
  values <- tolower(.clean_tax_value(tax[[domain_col[[1L]]]]))
  values <- values[!is.na(values)]
  if (!length(values)) return(FALSE)
  mean(values %in% c("bacteria", "archaea", "k__bacteria", "k__archaea", "d__bacteria", "d__archaea")) >= 0.5
}

.consensus_taxonomy <- function(tax, groups, features) {
  ranks <- .rank_columns(tax)
  if (!length(ranks)) return(NULL)
  out <- as.data.frame(matrix(NA_character_, nrow = length(features), ncol = length(ranks)),
                       stringsAsFactors = FALSE)
  names(out) <- names(ranks)
  rownames(out) <- features
  for (i in seq_along(features)) {
    rows <- which(groups == features[[i]])
    if (!length(rows)) next
    for (rank in names(ranks)) {
      values <- unique(stats::na.omit(.clean_tax_value(tax[[ranks[[rank]]]][rows])))
      if (length(values) == 1L) out[i, rank] <- values
    }
  }
  out
}

.match_taxonomy_to_counts <- function(counts, tax, allow_partial = FALSE) {
  if (!is.data.frame(tax) && !is.matrix(tax)) return(NULL)
  tax <- as.data.frame(tax, stringsAsFactors = FALSE, check.names = FALSE)
  features <- colnames(counts)
  if (is.null(features)) return(NULL)

  rn <- rownames(tax)
  if (!is.null(rn) && all(features %in% rn)) {
    return(tax[match(features, rn), , drop = FALSE])
  }

  taxa_col <- names(tax)[tolower(names(tax)) == "taxa"]
  if (length(taxa_col)) {
    groups <- as.character(tax[[taxa_col[[1L]]]])
    matched <- match(features, groups)
    if (all(!is.na(matched)) || (isTRUE(allow_partial) && any(!is.na(matched)))) {
      if (!anyDuplicated(groups)) {
        aligned <- tax[matched, , drop = FALSE]
        rownames(aligned) <- features
        return(aligned)
      }
      return(.consensus_taxonomy(tax, groups, features))
    }
  }
  NULL
}

.prepare_picrust_taxonomy <- function(counts, tax) {
  features <- colnames(counts)
  if (is.data.frame(tax) && "original_feature_id" %in% names(tax) &&
      all(features %in% tax$original_feature_id)) {
    aligned <- tax[match(features, tax$original_feature_id), , drop = FALSE]
  } else {
    aligned <- .match_taxonomy_to_counts(counts, tax)
  }
  if (is.null(aligned)) {
    stop(
      "PICRUSt2 requires taxonomy identifiers that match every ASV in the count table.",
      call. = FALSE
    )
  }
  aligned <- as.data.frame(aligned, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(aligned) != ncol(counts)) {
    stop("PICRUSt2 taxonomy and count tables have different ASV counts.", call. = FALSE)
  }
  rownames(aligned) <- colnames(counts)
  aligned$original_feature_id <- NULL
  data.frame(
    original_feature_id = colnames(counts),
    aligned,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.taxonomy_lineages <- function(tax) {
  tax <- as.data.frame(tax, stringsAsFactors = FALSE, check.names = FALSE)
  ranks <- .rank_columns(tax)
  if (!length(ranks)) {
    taxa_col <- names(tax)[tolower(names(tax)) == "taxa"]
    if (!length(taxa_col)) return(rep(NA_character_, nrow(tax)))
    values <- trimws(as.character(tax[[taxa_col[[1L]]]]))
    values[!nzchar(values) | is.na(values)] <- NA_character_
    full <- grepl("(^|;)\\s*[dkpcofgst]__", values, ignore.case = TRUE)
    prefixed <- grepl("^(uc_)?[kpcofgs]_", values, ignore.case = TRUE)
    values[prefixed] <- sub("^(?:uc_)?([kpcofgs])_", "\\1__", values[prefixed],
                            ignore.case = TRUE, perl = TRUE)
    values[!full & !prefixed] <- NA_character_
    return(values)
  }
  prefixes <- c(Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__",
                Family = "f__", Genus = "g__", Species = "s__", Strain = "t__")
  out <- character(nrow(tax))
  for (i in seq_len(nrow(tax))) {
    parts <- character()
    for (rank in names(ranks)) {
      value <- .clean_tax_value(tax[[ranks[[rank]]]][[i]])
      if (!is.na(value) && rank %in% c("Species", "Strain")) {
        value <- gsub("_", " ", value, fixed = TRUE)
      }
      if (!is.na(value)) parts <- c(parts, paste0(prefixes[[rank]], value))
    }
    out[[i]] <- paste(parts, collapse = ";")
  }
  out[!nzchar(out)] <- NA_character_
  out
}

.discover_picrust_sources <- function(study_dir) {
  counts_files <- list.files(
    study_dir,
    pattern = "_dada2_counts\\.rds\\.zip$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  parse_file <- file.path(study_dir, "parse.R")
  if (!file.exists(parse_file)) {
    study <- basename(normalizePath(study_dir, mustWork = TRUE))
    bundled <- tryCatch(.mutt_parser_file(study), error = function(e) "")
    if (nzchar(bundled)) parse_file <- bundled
  }
  if (file.exists(parse_file)) {
    parser_lines <- readLines(parse_file, warn = FALSE)
    parser_lines <- sub("#.*$", "", parser_lines)
    active_parser <- paste(parser_lines, collapse = "\n")
    counts_files <- counts_files[vapply(
      counts_files,
      function(path) grepl(basename(path), active_parser, fixed = TRUE),
      logical(1)
    )]
  }
  out <- lapply(sort(counts_files), function(counts_file) {
    prefix <- sub("_dada2_counts\\.rds\\.zip$", "", basename(counts_file), ignore.case = TRUE)
    tax_file <- file.path(study_dir, paste0(prefix, "_dada2_taxa.rds.zip"))
    list(
      id = prefix,
      type = "dada2",
      modality = "amplicon",
      counts_file = counts_file,
      tax_file = if (file.exists(tax_file)) tax_file else ""
    )
  })

  explicit_directories <- c(file.path(study_dir, "orientation_repair"), study_dir)
  explicit_required <- c(
    "asv_counts_by_run.tsv", "asv_sequences.tsv", "asv_orientation_qc.tsv"
  )
  explicit_directories <- explicit_directories[
    vapply(
      explicit_directories,
      function(path) all(vapply(
        file.path(path, explicit_required),
        function(candidate) nzchar(.resolve_study_file(candidate, must_work = FALSE)),
        logical(1)
      )),
      logical(1)
    )
  ]
  if (length(explicit_directories)) {
    explicit_directory <- explicit_directories[[1L]]
    out[[length(out) + 1L]] <- list(
      id = "asv",
      type = "explicit",
      modality = "amplicon",
      counts_file = .resolve_study_file(file.path(explicit_directory, "asv_counts_by_run.tsv")),
      sequences_file = .resolve_study_file(file.path(explicit_directory, "asv_sequences.tsv")),
      tax_file = .resolve_study_file(
        file.path(explicit_directory, "taxonomy_rdp19_toGenus.tsv"),
        must_work = FALSE
      ),
      orientation_qc_file = .resolve_study_file(
        file.path(explicit_directory, "asv_orientation_qc.tsv")
      )
    )
  }
  out
}

.load_picrust_source <- function(source) {
  if (source$type == "dada2") {
    counts <- .read_rds_zip(source$counts_file, "_counts\\.rds$")
    counts <- .normalize_count_matrix(counts, require_sequences = TRUE)
    sequences <- setNames(colnames(counts), colnames(counts))
    tax <- NULL
    if (nzchar(source$tax_file)) {
      tax <- as.data.frame(.read_rds_zip(source$tax_file, "_taxa\\.rds$"), stringsAsFactors = FALSE)
      if (is.null(rownames(tax)) && nrow(tax) == ncol(counts)) rownames(tax) <- colnames(counts)
    }
    repaired <- .collapse_reverse_complement_duplicates(counts, sequences, tax)
    counts <- repaired$counts
    sequences <- repaired$sequences
    tax <- repaired$taxonomy
  } else {
    raw_counts <- .read_study_delim(
      source$counts_file, check.names = FALSE, stringsAsFactors = FALSE
    )
    if (!"SampleID" %in% names(raw_counts)) stop("Explicit ASV table requires `SampleID`.", call. = FALSE)
    sample_ids <- raw_counts$SampleID
    raw_counts$SampleID <- NULL
    counts <- as.matrix(raw_counts)
    rownames(counts) <- sample_ids

    seq_map <- .read_study_delim(
      source$sequences_file, check.names = FALSE, stringsAsFactors = FALSE
    )
    if (!all(c("ASV", "Sequence") %in% names(seq_map))) {
      stop("Explicit sequence table requires `ASV` and `Sequence`.", call. = FALSE)
    }
    if (any(is.na(seq_map$ASV) | !nzchar(seq_map$ASV)) || anyDuplicated(seq_map$ASV)) {
      stop("Explicit sequence table requires unique, non-empty ASV IDs.", call. = FALSE)
    }
    if (any(is.na(seq_map$Sequence) | !nzchar(seq_map$Sequence)) ||
        anyDuplicated(toupper(seq_map$Sequence))) {
      stop("Explicit sequence table requires unique, non-empty sequence strings.", call. = FALSE)
    }
    sequences <- setNames(as.character(seq_map$Sequence), as.character(seq_map$ASV))
    if (!setequal(colnames(counts), names(sequences))) {
      stop("ASV counts and sequences do not match.", call. = FALSE)
    }
    sequences <- sequences[colnames(counts)]
    counts <- .normalize_count_matrix(counts)
    sequences <- sequences[colnames(counts)]
    if (any(is.na(sequences)) || .dna_fraction(sequences) < 1) {
      stop("Explicit ASV sequences contain missing or invalid DNA strings.", call. = FALSE)
    }
    if (!is.null(source$orientation_qc_file)) {
      if (!file.exists(source$orientation_qc_file)) {
        stop(
          "Explicit Caenepeel ASV input is missing `asv_orientation_qc.tsv`; ",
          "rebuild it with the orientation-normalized DADA2 workflow.",
          call. = FALSE
        )
      }
      orientation_qc <- .read_study_delim(
        source$orientation_qc_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      required_qc <- c(
        "asv_count", "exact_reverse_complement_pair_count", "all_asvs_forward_oriented"
      )
      if (nrow(orientation_qc) != 1L || !all(required_qc %in% names(orientation_qc)) ||
          orientation_qc$exact_reverse_complement_pair_count[[1L]] != 0 ||
          !isTRUE(orientation_qc$all_asvs_forward_oriented[[1L]]) ||
          orientation_qc$asv_count[[1L]] != length(sequences)) {
        stop("Caenepeel orientation QC did not certify the ASV input.", call. = FALSE)
      }
    }

    tax <- NULL
    if (file.exists(source$tax_file)) {
      tax <- .read_study_delim(
        source$tax_file, check.names = FALSE, stringsAsFactors = FALSE
      )
      id_col <- names(tax)[tolower(names(tax)) == "asv"]
      if (length(id_col)) rownames(tax) <- tax[[id_col[[1L]]]]
    }
  }

  .assert_no_reverse_complement_duplicates(sequences, "PICRUSt2 sequence input")

  aligned_taxonomy <- .prepare_picrust_taxonomy(counts, tax)

  list(
    counts = counts,
    sequences = sequences,
    taxonomy = aligned_taxonomy,
    prokaryotic = .is_prokaryotic_taxonomy(aligned_taxonomy),
    input_hash = .files_md5(c(
      source$counts_file,
      if (!is.null(source$sequences_file)) source$sequences_file,
      if (!is.null(source$orientation_qc_file)) source$orientation_qc_file,
      if (!is.null(source$tax_file) && file.exists(source$tax_file)) source$tax_file
    ))
  )
}

.read_stockholm_queries <- function(path, query_ids) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^(#|//|[[:space:]]*$)", lines)]
  fields <- strsplit(trimws(lines), "[[:space:]]+")
  keep <- lengths(fields) >= 2L
  ids <- vapply(fields[keep], `[[`, character(1), 1L)
  aligned <- vapply(fields[keep], `[[`, character(1), 2L)
  query <- ids %in% query_ids
  if (!any(query)) stop("HMM alignment did not return any query sequences.", call. = FALSE)
  combined <- vapply(
    split(aligned[query], ids[query]),
    paste0,
    collapse = "",
    FUN.VALUE = character(1)
  )
  missing <- setdiff(query_ids, names(combined))
  if (length(missing)) {
    stop("HMM alignment omitted ", length(missing), " query sequence(s).", call. = FALSE)
  }
  combined[query_ids]
}

.hmm_orientation_fractions <- function(sequences, reverse_sequences, reference_msa,
                                       reference_hmm, hmmalign, prefix, directory) {
  n <- length(sequences)
  forward_ids <- paste0("forward_", seq_len(n))
  reverse_ids <- paste0("reverse_", seq_len(n))
  ids <- c(forward_ids, reverse_ids)
  query_sequences <- c(unname(sequences), unname(reverse_sequences))
  fasta <- file.path(directory, paste0(prefix, "queries.fna"))
  stockholm <- file.path(directory, paste0(prefix, "alignment.stockholm"))
  stdout <- file.path(directory, paste0(prefix, "hmmalign.stdout.log"))
  stderr <- file.path(directory, paste0(prefix, "hmmalign.stderr.log"))
  writeLines(as.vector(rbind(paste0(">", ids), query_sequences)), fasta)
  run <- .run_system2(
    hmmalign,
    c(
      "--trim", "--dna", "--mapali", shQuote(reference_msa),
      "--informat", "FASTA", "-o", shQuote(stockholm),
      shQuote(reference_hmm), shQuote(fasta)
    ),
    stdout,
    stderr
  )
  if (run$status != 0L || !file.exists(stockholm)) {
    stop("PICRUSt2 orientation HMM alignment failed; see ", stderr, call. = FALSE)
  }
  aligned <- .read_stockholm_queries(stockholm, ids)
  aligned_length <- nchar(gsub("[-.]", "", aligned))
  fractions <- as.numeric(aligned_length / nchar(query_sequences))
  list(forward = fractions[seq_len(n)], reverse = fractions[n + seq_len(n)])
}

.orient_picrust_sequences <- function(sequences, directory, tools) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  dependencies <- tools$picrust2_dependencies
  hmmalign <- if (!is.null(dependencies) && "hmmalign" %in% names(dependencies)) {
    dependencies[["hmmalign"]]
  } else {
    ""
  }
  references <- tools$picrust2_references
  required_references <- c("bacteria_msa", "bacteria_hmm", "archaea_msa", "archaea_hmm")
  if (!nzchar(hmmalign)) {
    stop("PICRUSt2 orientation checking requires hmmalign on PATH.", call. = FALSE)
  }
  if (is.null(references) || !all(required_references %in% names(references)) ||
      any(!nzchar(references[required_references])) ||
      any(!file.exists(references[required_references]))) {
    stop("PICRUSt2 bacterial and archaeal reference HMM files were not found.", call. = FALSE)
  }

  reverse_sequences <- .reverse_complement_dna(sequences)
  bacteria <- .hmm_orientation_fractions(
    sequences, reverse_sequences,
    references[["bacteria_msa"]], references[["bacteria_hmm"]],
    hmmalign, "bacteria_", directory
  )
  archaea <- .hmm_orientation_fractions(
    sequences, reverse_sequences,
    references[["archaea_msa"]], references[["archaea_hmm"]],
    hmmalign, "archaea_", directory
  )

  forward_fraction <- pmax(bacteria$forward, archaea$forward)
  reverse_fraction <- pmax(bacteria$reverse, archaea$reverse)
  use_reverse <- reverse_fraction > forward_fraction
  oriented <- sequences
  oriented[use_reverse] <- reverse_sequences[use_reverse]
  chosen_fraction <- ifelse(use_reverse, reverse_fraction, forward_fraction)
  audit <- data.frame(
    original_feature_id = names(sequences),
    orientation = ifelse(use_reverse, "reverse_complement", "supplied"),
    forward_aligned_fraction = forward_fraction,
    reverse_aligned_fraction = reverse_fraction,
    chosen_aligned_fraction = chosen_fraction,
    passes_min_align = chosen_fraction >= .picrust_min_align,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    audit,
    file.path(directory, "asv_orientation.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  list(sequences = oriented, audit = audit)
}

.write_picrust_inputs <- function(counts, sequences, directory, biom_command) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (!nzchar(biom_command)) {
    stop("The BIOM command-line tool is required to prepare PICRUSt2 input.", call. = FALSE)
  }
  ids <- paste0("ASV", seq_len(ncol(counts)))
  seq_values <- if (!is.null(names(sequences)) && all(colnames(counts) %in% names(sequences))) {
    unname(sequences[colnames(counts)])
  } else {
    colnames(counts)
  }
  if (any(is.na(seq_values)) || .dna_fraction(seq_values) < 1) {
    stop("PICRUSt2 sequences are missing or invalid.", call. = FALSE)
  }

  fasta <- file.path(directory, "study_sequences.fna")
  writeLines(as.vector(rbind(paste0(">", ids), toupper(seq_values))), fasta)

  mapping <- data.frame(
    picrust_id = ids,
    original_feature_id = colnames(counts),
    sequence = toupper(seq_values),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  mapping_path <- file.path(directory, "asv_id_map.tsv")
  utils::write.table(
    mapping, mapping_path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
  )

  table <- t(counts)
  rownames(table) <- ids
  classic_path <- file.path(directory, "asv_counts.tsv")
  utils::write.table(
    cbind(`#OTU ID` = rownames(table), as.data.frame(table, check.names = FALSE)),
    classic_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  biom_path <- file.path(directory, "asv_counts.biom")
  convert <- .run_system2(
    biom_command,
    c("convert", "-i", shQuote(classic_path), "-o", shQuote(biom_path),
      "--table-type", shQuote("OTU table"), "--to-hdf5"),
    file.path(directory, "biom_convert.stdout.log"),
    file.path(directory, "biom_convert.stderr.log")
  )
  if (convert$status != 0L || !file.exists(biom_path)) {
    stop("Failed to convert the PICRUSt2 count table to BIOM format.", call. = FALSE)
  }
  validate <- .run_system2(
    biom_command,
    c("validate-table", "-i", shQuote(biom_path)),
    file.path(directory, "biom_validate.stdout.log"),
    file.path(directory, "biom_validate.stderr.log")
  )
  if (validate$status != 0L) {
    stop("Generated PICRUSt2 BIOM input did not pass `biom validate-table`.", call. = FALSE)
  }

  list(
    fasta = fasta,
    table = biom_path,
    classic_table = classic_path,
    mapping = mapping_path,
    mapping_data = mapping
  )
}

.read_function_table <- function(
  path,
  expected_samples = NULL,
  allow_missing_samples = FALSE
) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  if (!length(lines)) stop("Functional output is empty: ", path, call. = FALSE)

  tab_lines <- which(grepl("\t", lines, fixed = TRUE))
  if (!length(tab_lines)) stop("Functional output is not tab-delimited: ", path, call. = FALSE)
  header <- tab_lines[[1L]]
  if (!is.null(expected_samples)) {
    expected_headers <- vapply(tab_lines, function(i) {
      fields <- strsplit(lines[[i]], "\t", fixed = TRUE)[[1L]]
      all(expected_samples %in% fields)
    }, logical(1))
    if (any(expected_headers)) header <- tab_lines[which(expected_headers)[[1L]]]
  }
  text <- paste(lines[header:length(lines)], collapse = "\n")
  text_connection <- textConnection(text)
  on.exit(close(text_connection), add = TRUE)
  x <- utils::read.delim(text_connection, check.names = FALSE, stringsAsFactors = FALSE,
                         comment.char = "", quote = "")
  if (ncol(x) < 2L) stop("Functional output has no sample columns: ", path, call. = FALSE)
  feature_ids <- as.character(x[[1L]])
  x[[1L]] <- NULL
  numeric_cols <- vapply(x, function(z) all(!is.na(suppressWarnings(as.numeric(z)))), logical(1))
  x <- x[, numeric_cols, drop = FALSE]
  if (!ncol(x)) stop("Functional output has no numeric sample columns: ", path, call. = FALSE)
  x[] <- lapply(x, as.numeric)
  out <- t(as.matrix(x))
  colnames(out) <- feature_ids
  if (!is.null(expected_samples)) {
    expected_samples <- as.character(expected_samples)
    observed_samples <- rownames(out)
    if (anyNA(observed_samples) || any(!nzchar(observed_samples)) ||
        anyDuplicated(observed_samples)) {
      stop("Functional output requires unique, non-empty sample IDs.", call. = FALSE)
    }
    missing_samples <- setdiff(expected_samples, observed_samples)
    extra_samples <- setdiff(observed_samples, expected_samples)
    if (length(extra_samples) || (length(missing_samples) && !allow_missing_samples)) {
      stop(
        "Functional output sample IDs do not match its input counts: ",
        length(missing_samples), " expected sample(s) missing and ",
        length(extra_samples), " unexpected sample(s) present.",
        call. = FALSE
      )
    }
    if (length(missing_samples)) {
      restored <- matrix(
        0,
        nrow = length(expected_samples),
        ncol = ncol(out),
        dimnames = list(expected_samples, colnames(out))
      )
      restored[observed_samples, ] <- out
      out <- restored
    } else {
      out <- out[expected_samples, , drop = FALSE]
    }
    if (isTRUE(allow_missing_samples)) {
      attr(out, "sample_reconciliation") <- list(
        expected_samples = expected_samples,
        observed_samples = observed_samples,
        missing_samples = missing_samples,
        extra_samples = extra_samples,
        zero_filled_missing_samples = length(missing_samples)
      )
    }
  }
  out
}

.align_picrust_pathway_samples <- function(pathway_table, expected_samples) {
  expected_samples <- as.character(expected_samples)
  observed_samples <- rownames(pathway_table)
  missing_samples <- setdiff(expected_samples, observed_samples)
  extra_samples <- setdiff(observed_samples, expected_samples)

  if (length(extra_samples)) {
    stop(
      "PICRUSt2 pathway output contains ", length(extra_samples),
      " unexpected sample(s): ", paste(extra_samples, collapse = ", "),
      call. = FALSE
    )
  }

  warning_message <- ""
  if (length(missing_samples)) {
    zero_rows <- matrix(
      0,
      nrow = length(missing_samples),
      ncol = ncol(pathway_table),
      dimnames = list(missing_samples, colnames(pathway_table))
    )
    pathway_table <- rbind(pathway_table, zero_rows)
    warning_message <- sprintf(
      paste0(
        "Pathway output omitted %d of %d samples; ",
        "missing sample%s %s zero-filled."
      ),
      length(missing_samples),
      length(expected_samples),
      if (length(missing_samples) == 1L) "" else "s",
      if (length(missing_samples) == 1L) "was" else "were"
    )
  }

  pathway_table <- pathway_table[expected_samples, , drop = FALSE]
  list(
    table = pathway_table,
    warning = warning_message,
    zero_filled_samples = missing_samples
  )
}

.find_one <- function(directory, pattern, required = TRUE) {
  hits <- list.files(directory, pattern = pattern, recursive = TRUE, full.names = TRUE,
                     ignore.case = TRUE)
  if (!length(hits)) {
    if (required) stop("Expected functional output not found: ", pattern, call. = FALSE)
    return("")
  }
  hits[[1L]]
}

.picrust_contribution_descriptor <- function(raw_dir, path) {
  if (!length(path) || is.na(path) || !file.exists(path)) {
    stop("Expected PICRUSt2 contribution output was not found.", call. = FALSE)
  }
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header)) stop("PICRUSt2 contribution output is empty: ", path, call. = FALSE)
  columns <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  if (!all(c("sample", "function", "taxon") %in% columns)) {
    stop("PICRUSt2 contribution output has an unexpected schema: ", path, call. = FALSE)
  }
  list(
    relative_path = file.path("raw", substring(path, nchar(raw_dir) + 2L)),
    columns = columns,
    size_bytes = unname(file.info(path)$size)
  )
}

.parse_picrust_output <- function(raw_dir, counts) {
  all_pred <- list.files(raw_dir, pattern = "^pred_metagenome_unstrat\\.tsv\\.gz$",
                         recursive = TRUE, full.names = TRUE)
  ec_path <- all_pred[grepl("EC_metagenome_out", all_pred, fixed = TRUE)][1L]
  ko_path <- all_pred[grepl("KO_metagenome_out", all_pred, fixed = TRUE)][1L]
  if (!length(ec_path) || is.na(ec_path)) stop("PICRUSt2 EC output is missing.", call. = FALSE)
  if (!length(ko_path) || is.na(ko_path)) stop("PICRUSt2 KO output is missing.", call. = FALSE)
  abundance_path <- .find_one(raw_dir, "^path_abun_unstrat\\.tsv\\.gz$")
  coverage_path <- .find_one(raw_dir, "^path_cov_unstrat\\.tsv\\.gz$")
  nsti_path <- .find_one(
    raw_dir,
    "(weighted_nsti|marker.*nsti|nsti.*predicted).*\\.tsv(\\.gz)?$",
    required = FALSE
  )

  nsti <- NULL
  if (nzchar(nsti_path)) {
    con <- if (grepl("\\.gz$", nsti_path)) gzfile(nsti_path, "rt") else file(nsti_path, "rt")
    nsti <- utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
    close(con)
  }
  expected <- rownames(counts)
  ec <- .read_function_table(ec_path, expected, allow_missing_samples = TRUE)
  ko <- .read_function_table(ko_path, expected, allow_missing_samples = TRUE)
  metacyc_abundance <- .read_function_table(abundance_path)
  metacyc_coverage <- .read_function_table(coverage_path)
  ec_reconciliation <- attr(ec, "sample_reconciliation", exact = TRUE)
  ko_reconciliation <- attr(ko, "sample_reconciliation", exact = TRUE)
  ec_samples <- ec_reconciliation$observed_samples
  ko_samples <- ko_reconciliation$observed_samples
  if (!identical(ec_samples, ko_samples)) {
    stop(
      "PICRUSt2 EC and KO output tables disagree about retained sample IDs.",
      call. = FALSE
    )
  }
  abundance_alignment <- .align_picrust_pathway_samples(
    metacyc_abundance, ec_samples
  )
  coverage_alignment <- .align_picrust_pathway_samples(
    metacyc_coverage, ec_samples
  )
  pathway_zero_filled_samples <- unique(c(
    abundance_alignment$zero_filled_samples,
    coverage_alignment$zero_filled_samples
  ))
  warning_messages <- unique(Filter(
    nzchar,
    c(abundance_alignment$warning, coverage_alignment$warning)
  ))
  picrust_warning <- paste(warning_messages, collapse = " ")
  metacyc_abundance <- abundance_alignment$table
  metacyc_coverage <- coverage_alignment$table
  missing_samples <- ec_reconciliation$missing_samples
  if (length(missing_samples)) {
    restore <- function(x) {
      out <- matrix(0, nrow = length(expected), ncol = ncol(x),
                    dimnames = list(expected, colnames(x)))
      out[ec_samples, ] <- x
      out
    }
    metacyc_abundance <- restore(metacyc_abundance)
    metacyc_coverage <- restore(metacyc_coverage)
  }
  attr(ec, "sample_reconciliation") <- NULL
  attr(ko, "sample_reconciliation") <- NULL
  contribution_paths <- list.files(
    raw_dir,
    pattern = "^(pred_metagenome|path_abun)_contrib\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  ec_contribution <- contribution_paths[
    grepl("EC_metagenome_out", contribution_paths, fixed = TRUE)
  ][1L]
  ko_contribution <- contribution_paths[
    grepl("KO_metagenome_out", contribution_paths, fixed = TRUE)
  ][1L]
  pathway_contribution <- contribution_paths[
    grepl("path_abun_contrib.tsv.gz", contribution_paths, fixed = TRUE)
  ][1L]
  stratified <- list(
    ec = .picrust_contribution_descriptor(raw_dir, ec_contribution),
    ko = .picrust_contribution_descriptor(raw_dir, ko_contribution),
    metacyc_abundance = .picrust_contribution_descriptor(raw_dir, pathway_contribution)
  )

  list(
    ec = ec,
    ko = ko,
    metacyc_abundance = metacyc_abundance,
    metacyc_coverage = metacyc_coverage,
    stratified = stratified,
    nsti = nsti,
    qc = list(
      input_samples = nrow(counts),
      input_asvs = ncol(counts),
      ec_features = ncol(ec),
      ko_features = ncol(ko),
      metacyc_pathways = ncol(metacyc_abundance),
      nsti_records = if (is.null(nsti)) NA_integer_ else nrow(nsti),
      picrust_output_samples = length(ec_samples),
      zero_filled_missing_samples = length(missing_samples),
      pathway_zero_filled_samples = length(pathway_zero_filled_samples)
    ),
    sample_reconciliation = data.frame(
      sample_id = expected,
      status = ifelse(expected %in% missing_samples, "zero_filled_missing", "retained"),
      stringsAsFactors = FALSE
    ),
    picrust_warning = picrust_warning,
    pathway_zero_filled_samples = pathway_zero_filled_samples
  )
}

#' Read stratified PICRUSt2 ASV contributions
#'
#' Reads one file-backed PICRUSt2 contribution table from a branch returned by
#' `mutt(functional = TRUE)`. PICRUSt2's internal sequence identifiers are
#' retained in `picrust_id`, and `original_feature_id` identifies the ASV in the
#' source count and taxonomy tables. Taxonomy columns are joined to every
#' contribution row. The full table is read by default and can be large.
#'
#' @param x A single PICRUSt2 branch result, for example
#'   `repo[[study]][["function"]]$picrust2[[branch]]`.
#' @param type One of `"ec"`, `"ko"`, or `"metacyc_abundance"`.
#' @param n_max Maximum number of contribution rows to read. Use a small value
#'   for inspection; the default reads all rows.
#'
#' @return A data frame with sample, function, original ASV identifier, and
#'   PICRUSt2 contribution quantities.
#' @noRd
read_picrust2_contributions <- function(
  x,
  type = c("ec", "ko", "metacyc_abundance"),
  n_max = Inf
) {
  type <- match.arg(type)
  if (!is.list(x) || is.null(x$stratified) || is.null(x$asv_mapping) ||
      is.null(x$provenance$output_directory)) {
    stop("`x` must be a PICRUSt2 branch result returned by MUTT.", call. = FALSE)
  }
  if (!is.numeric(n_max) || length(n_max) != 1L || is.na(n_max) || n_max <= 0) {
    stop("`n_max` must be a positive number or Inf.", call. = FALSE)
  }
  descriptor <- x$stratified[[type]]
  path <- file.path(x$provenance$output_directory, descriptor$relative_path)
  if (!file.exists(path)) {
    stop(
      "The stratified PICRUSt2 file was not found at ", path,
      ". Re-run the study with functional = \"REBUILD\" if its cache moved.",
      call. = FALSE
    )
  }
  limit <- if (is.infinite(n_max)) Inf else as.integer(n_max)
  out <- readr::read_tsv(
    path,
    n_max = limit,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  taxon <- match("taxon", names(out))
  if (is.na(taxon)) stop("The contribution table has no `taxon` column.", call. = FALSE)
  names(out)[taxon] <- "picrust_id"
  mapped <- x$asv_mapping$original_feature_id[
    match(as.character(out$picrust_id), x$asv_mapping$picrust_id)
  ]
  if (anyNA(mapped)) {
    stop("A PICRUSt2 contribution identifier could not be mapped to its source ASV.", call. = FALSE)
  }
  out <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  out$original_feature_id <- mapped
  taxonomy_rows <- match(mapped, x$taxonomy$original_feature_id)
  if (anyNA(taxonomy_rows)) {
    stop("A contribution ASV could not be matched to the cached taxonomy table.", call. = FALSE)
  }
  taxonomy <- x$taxonomy[
    taxonomy_rows,
    setdiff(names(x$taxonomy), "original_feature_id"),
    drop = FALSE
  ]
  collisions <- names(taxonomy) %in% names(out)
  names(taxonomy)[collisions] <- paste0("taxonomy_", names(taxonomy)[collisions])
  out <- cbind(out, taxonomy)
  out
}

.run_system2 <- function(command, args, stdout, stderr) {
  started <- proc.time()[["elapsed"]]
  invocation <- .command_invocation(command, args)
  status <- tryCatch(
    system2(invocation$command, args = invocation$args, stdout = stdout, stderr = stderr),
    error = function(e) structure(127L, message = conditionMessage(e))
  )
  list(status = as.integer(status), runtime = proc.time()[["elapsed"]] - started)
}

.cache_metadata <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

.cache_is_valid <- function(directory, input_hash, tool_version = "", database = "") {
  metadata <- .cache_metadata(file.path(directory, "manifest.json"))
  if (!is.list(metadata) || !identical(as.character(metadata$input_hash), input_hash) ||
      !identical(as.character(metadata$engine_version), .functional_engine_version) ||
      !file.exists(file.path(directory, "result.rds"))) {
    return(FALSE)
  }
  if (nzchar(tool_version) && !identical(as.character(metadata$tool_version), tool_version)) {
    return(FALSE)
  }
  if (nzchar(database) && !identical(as.character(metadata$database), database)) {
    return(FALSE)
  }
  TRUE
}

.restore_functional_archive <- function(study_dir) {
  target <- file.path(study_dir, "functional")
  archive <- paste0(target, ".zip")
  if (dir.exists(target) || !file.exists(archive)) return(FALSE)

  members <- utils::unzip(archive, list = TRUE)$Name
  normalized <- gsub("\\\\", "/", members)
  components <- strsplit(normalized, "/", fixed = TRUE)
  unsafe <- !length(members) ||
    any(grepl("^/|^[A-Za-z]:", normalized)) ||
    any(vapply(components, function(x) any(x == ".."), logical(1))) ||
    any(!(normalized == "functional" | grepl("^functional/", normalized)))
  if (unsafe) {
    stop("Functional cache archive contains an unsafe or unexpected path: ", archive,
         call. = FALSE)
  }

  staging <- tempfile("functional-restore-", tmpdir = study_dir)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(archive, exdir = staging)
  restored <- file.path(staging, "functional")
  if (!dir.exists(restored)) {
    stop("Functional cache archive did not contain a `functional/` directory: ", archive,
         call. = FALSE)
  }
  restored_files <- list.files(restored, recursive = TRUE, full.names = TRUE,
                               all.files = TRUE, no.. = TRUE)
  if (length(restored_files) && any(nzchar(Sys.readlink(restored_files)))) {
    stop("Functional cache archives may not contain symbolic links: ", archive,
         call. = FALSE)
  }
  if (!file.rename(restored, target)) {
    stop("Could not publish restored functional cache: ", archive, call. = FALSE)
  }
  TRUE
}

.compact_picrust_cache <- function(directory) {
  removable <- c(
    file.path(directory, "input", "orientation"),
    file.path(directory, "raw", "intermediate")
  )
  for (path in removable[dir.exists(removable)]) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  invisible(directory)
}

.retain_failed_picrust_staging <- function(staging, target) {
  failed <- paste0(target, ".failed")
  if (dir.exists(failed) || file.exists(failed)) {
    failed <- paste0(
      failed, "-", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
      "-", Sys.getpid()
    )
  }
  dir.create(dirname(failed), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging, failed)) return(staging)
  failed
}

.publish_cache <- function(staging, target) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  backup <- ""
  if (dir.exists(target)) {
    backup <- paste0(target, ".previous")
    if (dir.exists(backup)) unlink(backup, recursive = TRUE)
    if (!file.rename(target, backup)) stop("Could not protect the previous functional result.", call. = FALSE)
  }
  if (!file.rename(staging, target)) {
    if (nzchar(backup) && dir.exists(backup)) file.rename(backup, target)
    stop("Could not publish the functional result cache.", call. = FALSE)
  }
  if (nzchar(backup) && dir.exists(backup)) unlink(backup, recursive = TRUE)
  invisible(target)
}

.run_picrust2 <- function(counts, sequences, branch, source, functional_dir, mode, tools,
                          input_hash = NULL, taxonomy = NULL) {
  target <- file.path(functional_dir, "picrust2", .safe_branch_name(branch))
  taxonomy <- .prepare_picrust_taxonomy(counts, taxonomy)
  if (is.null(input_hash)) {
    input_hash <- .object_md5(list(counts = counts, sequences = sequences, taxonomy = taxonomy))
  }
  if (mode != "rebuild" && .cache_is_valid(target, input_hash, tools$picrust2_version)) {
    result <- readRDS(file.path(target, "result.rds"))
    class(result) <- unique(c("mutt_picrust_branch", class(result)))
    result$provenance$output_directory <- target
    metadata <- .cache_metadata(file.path(target, "manifest.json"))
    return(list(
      result = result,
      manifest = .manifest_row("picrust2", branch, "cached", source = source,
                               samples = nrow(counts), features = ncol(counts),
                               input_hash = input_hash, tool_version = metadata$tool_version,
                               runtime_seconds = 0, output_directory = target)
    ))
  }
  if (!nzchar(tools$picrust2)) {
    return(list(result = NULL, manifest = .manifest_row(
      "picrust2", branch, "skipped", "picrust2_pipeline.py was not found on PATH.", source,
      nrow(counts), ncol(counts), input_hash, output_directory = target
    )))
  }

  required <- c("place_seqs", "hsp", "metagenome", "pathway", "biom", "hmmalign")
  dependencies <- tools$picrust2_dependencies
  missing_dependencies <- if (is.null(dependencies)) {
    character()
  } else {
    required[!nzchar(dependencies[required]) | is.na(dependencies[required])]
  }
  if (length(missing_dependencies)) {
    stop(
      "PICRUSt2 is incomplete; missing required command(s): ",
      paste(missing_dependencies, collapse = ", "),
      call. = FALSE
    )
  }

  biom_command <- if (!is.null(dependencies) && nzchar(dependencies[["biom"]])) {
    dependencies[["biom"]]
  } else {
    unname(Sys.which("biom"))
  }
  processes <- .functional_processes()

  staging <- tempfile(paste0(.safe_branch_name(branch), "_"), tmpdir = functional_dir)
  dir.create(staging, recursive = TRUE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  orientation <- .orient_picrust_sequences(
    sequences,
    file.path(staging, "input", "orientation"),
    tools
  )
  if (!any(orientation$audit$passes_min_align)) {
    stop(
      "No ASV passed the PICRUSt2 ", .picrust_min_align,
      " alignment threshold in either orientation.",
      call. = FALSE
    )
  }
  inputs <- .write_picrust_inputs(
    counts, orientation$sequences, file.path(staging, "input"), biom_command = biom_command
  )
  orientation_rows <- match(
    inputs$mapping_data$original_feature_id,
    orientation$audit$original_feature_id
  )
  inputs$mapping_data <- cbind(
    inputs$mapping_data,
    orientation$audit[
      orientation_rows,
      setdiff(names(orientation$audit), "original_feature_id"),
      drop = FALSE
    ]
  )
  inputs$mapping_data$taxonomy_row_id <- inputs$mapping_data$original_feature_id
  utils::write.table(
    inputs$mapping_data,
    inputs$mapping,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  raw_dir <- file.path(staging, "raw")
  logs_dir <- file.path(staging, "logs")
  dir.create(logs_dir)
  old_path <- Sys.getenv("PATH", unset = "")
  tool_bin <- dirname(tools$picrust2)
  if (nzchar(tool_bin)) {
    withr::local_envvar(PATH = paste(tool_bin, old_path, sep = .Platform$path.sep))
  }
  run <- .run_system2(
    tools$picrust2,
    c(
      "-s", shQuote(inputs$fasta),
      "-i", shQuote(inputs$table),
      "-o", shQuote(raw_dir),
      "-p", as.character(processes),
      "--min_align", as.character(.picrust_min_align),
      "--coverage",
      "--stratified",
      "--verbose"
    ),
    file.path(logs_dir, "stdout.log"),
    file.path(logs_dir, "stderr.log")
  )
  if (run$status != 0L) {
    failed <- .retain_failed_picrust_staging(staging, target)
    log_location <- file.path(failed, "logs")
    stop("PICRUSt2 failed for ", branch, "; see ", log_location, call. = FALSE)
  }
  result <- tryCatch(
    .parse_picrust_output(raw_dir, counts),
    error = function(error) {
      failed <- .retain_failed_picrust_staging(staging, target)
      stop(
        "PICRUSt2 completed but output validation failed for ", branch,
        ": ", conditionMessage(error), " Raw outputs retained at ", failed, ".",
        call. = FALSE
      )
    }
  )
  result$asv_mapping <- inputs$mapping_data
  result$taxonomy <- taxonomy
  retained_ids <- if (is.null(result$nsti) || !nrow(result$nsti)) {
    character()
  } else {
    as.character(result$nsti[[1L]])
  }
  retained <- inputs$mapping_data$picrust_id %in% retained_ids
  input_reads <- sum(counts)
  retained_reads <- sum(counts[, retained, drop = FALSE])
  result$qc$retained_asvs <- sum(retained)
  result$qc$excluded_asvs <- sum(!retained)
  result$qc$retained_read_fraction <- if (input_reads > 0) retained_reads / input_reads else NA_real_
  result$qc$reverse_complemented_asvs <- sum(
    inputs$mapping_data$orientation == "reverse_complement"
  )
  result$qc$orientation_passed_asvs <- sum(inputs$mapping_data$passes_min_align)
  result$qc$orientation_failed_asvs <- sum(!inputs$mapping_data$passes_min_align)
  result$qc$min_align <- .picrust_min_align
  result$qc$processes <- processes
  if (nrow(result$sample_reconciliation)) {
    utils::write.table(
      result$sample_reconciliation,
      file.path(staging, "sample_reconciliation.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  if (length(result$pathway_zero_filled_samples)) {
    utils::write.table(
      data.frame(
        sample_id = result$pathway_zero_filled_samples,
        action = "zero_filled_in_pathway_output",
        stringsAsFactors = FALSE
      ),
      file.path(staging, "pathway_zero_filled_samples.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  result$provenance <- list(source = source, input_hash = input_hash,
                            tool_version = tools$picrust2_version,
                            biom_version = tools$biom_version,
                            output_directory = target,
                            orientation_method = "maximum aligned fraction across PICRUSt2 bacterial and archaeal reference HMMs",
                            min_align = .picrust_min_align,
                            processes = processes)
  class(result) <- unique(c("mutt_picrust_branch", class(result)))
  .compact_picrust_cache(staging)
  saveRDS(result, file.path(staging, "result.rds"))
  metadata <- list(method = "picrust2", branch = branch, source = source,
                   input_hash = input_hash, tool_version = tools$picrust2_version,
                   engine_version = .functional_engine_version,
                   biom_version = tools$biom_version,
                   processes = processes,
                   runtime_seconds = run$runtime)
  jsonlite::write_json(metadata, file.path(staging, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  .publish_cache(staging, target)
  picrust_status <- if (nzchar(result$picrust_warning)) {
    "generated_with_warning"
  } else {
    "generated"
  }
  picrust_reason <- if (nzchar(result$picrust_warning)) result$picrust_warning else ""
  list(result = result, manifest = .manifest_row(
    "picrust2", branch, picrust_status, reason = picrust_reason,
    source = source, samples = nrow(counts),
    features = ncol(counts), input_hash = input_hash,
    tool_version = tools$picrust2_version, runtime_seconds = run$runtime,
    output_directory = target
  ))
}

.write_faprotax_input <- function(counts, lineages, path) {
  table <- t(counts)
  out <- cbind(feature = rownames(table), as.data.frame(table, check.names = FALSE),
               taxonomy = unname(lineages))
  names(out)[[1L]] <- "#OTU ID"
  utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

.parse_assignment_report <- function(path) {
  if (!file.exists(path)) return(data.frame(line = character(), stringsAsFactors = FALSE))
  lines <- readLines(path, warn = FALSE)
  data.frame(line = lines, stringsAsFactors = FALSE)
}

.run_faprotax <- function(counts, taxonomy, branch, source, functional_dir, mode, tools) {
  lineages <- .taxonomy_lineages(taxonomy)
  keep <- !is.na(lineages) & colSums(counts) > 0
  counts <- counts[, keep, drop = FALSE]
  lineages <- lineages[keep]
  if (!ncol(counts)) {
    return(list(result = NULL, manifest = .manifest_row(
      "faprotax", branch, "skipped", "No classified, nonzero taxa were available.", source
    )))
  }

  target <- file.path(functional_dir, "faprotax", .safe_branch_name(branch))
  input_hash <- .object_md5(list(counts = counts, lineages = lineages))
  if (mode != "rebuild" && .cache_is_valid(
    target, input_hash, tools$faprotax_version, tools$faprotax_db_version
  )) {
    result <- readRDS(file.path(target, "result.rds"))
    result$provenance$database <- tools$faprotax_db
    metadata <- .cache_metadata(file.path(target, "manifest.json"))
    return(list(result = result, manifest = .manifest_row(
      "faprotax", branch, "cached", source = source, samples = nrow(counts),
      features = ncol(counts), input_hash = input_hash,
      tool_version = metadata$tool_version, runtime_seconds = 0,
      output_directory = target
    )))
  }
  if (!nzchar(tools$faprotax) || !nzchar(tools$faprotax_db)) {
    reason <- if (!nzchar(tools$faprotax)) {
      "collapse_table.py was not found on PATH."
    } else {
      "FAPROTAX.txt was not found beside the installed FAPROTAX executable."
    }
    return(list(result = NULL, manifest = .manifest_row(
      "faprotax", branch, "skipped", reason, source, nrow(counts), ncol(counts),
      input_hash, output_directory = target
    )))
  }

  staging <- tempfile(paste0(.safe_branch_name(branch), "_"), tmpdir = functional_dir)
  dir.create(staging, recursive = TRUE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  input_path <- file.path(staging, "taxon_counts.tsv")
  output_path <- file.path(staging, "functional_abundance.tsv")
  report_path <- file.path(staging, "assignments.txt")
  logs_dir <- file.path(staging, "logs")
  dir.create(logs_dir)
  .write_faprotax_input(counts, lineages, input_path)
  run <- .run_system2(
    tools$faprotax,
    c("-i", shQuote(input_path), "-o", shQuote(output_path),
      "-g", shQuote(tools$faprotax_db), "-d", "taxonomy", "-c", shQuote("#"),
      "--omit_columns", "0", "--column_names_are_in", "last_comment_line",
      "-r", shQuote(report_path), "-n", "columns_after_collapsing", "-v"),
    file.path(logs_dir, "stdout.log"),
    file.path(logs_dir, "stderr.log")
  )
  if (run$status != 0L) stop("FAPROTAX failed for ", branch, "; see ", logs_dir, call. = FALSE)
  abundance <- .read_function_table(output_path, rownames(counts))
  assignments <- .parse_assignment_report(report_path)
  result <- list(
    abundance = abundance,
    assignments = assignments,
    qc = list(
      input_samples = nrow(counts),
      input_taxa = ncol(counts),
      functional_groups = ncol(abundance)
    ),
    provenance = list(source = source, input_hash = input_hash,
                      tool_version = tools$faprotax_version,
                      database = tools$faprotax_db)
  )
  saveRDS(result, file.path(staging, "result.rds"))
  metadata <- list(method = "faprotax", branch = branch, source = source,
                   input_hash = input_hash, tool_version = tools$faprotax_version,
                   database = tools$faprotax_db_version,
                   database_path = tools$faprotax_db,
                   engine_version = .functional_engine_version,
                   runtime_seconds = run$runtime)
  jsonlite::write_json(metadata, file.path(staging, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  .publish_cache(staging, target)
  list(result = result, manifest = .manifest_row(
    "faprotax", branch, "generated", source = source, samples = nrow(counts),
    features = ncol(counts), input_hash = input_hash,
    tool_version = tools$faprotax_version, runtime_seconds = run$runtime,
    output_directory = target
  ))
}

.functional_faprotax_pairs <- function(parsed) {
  counts <- .flatten_tables(parsed$counts)
  proportions <- .flatten_tables(parsed$proportions)
  tax <- .flatten_tables(parsed$tax)
  out <- list()
  excluded <- .functional_faprotax_exclusion_template()
  for (path in names(tax)) {
    abundance <- NULL
    input_type <- ""
    if (path %in% names(counts)) {
      abundance <- tryCatch(.normalize_count_matrix(counts[[path]]), error = function(e) NULL)
      if (!is.null(abundance)) input_type <- "counts"
    }
    if (is.null(abundance) && path %in% names(proportions)) {
      abundance <- tryCatch(
        .normalize_abundance_matrix(proportions[[path]]),
        error = function(e) NULL
      )
      if (!is.null(abundance)) input_type <- "proportions"
    }
    if (is.null(abundance)) next
    source <- paste0(input_type, ":", path)
    modality <- .functional_branch_modality(path)
    if (!identical(modality, "amplicon")) {
      excluded <- rbind(
        excluded,
        data.frame(
          branch = path,
          reason = if (identical(modality, "metagenomic")) {
            "Excluded from FAPROTAX because this is a shotgun/metagenomic branch."
          } else {
            "Excluded from FAPROTAX because the branch is not explicitly identified as amplicon data."
          },
          source = source,
          samples = nrow(abundance),
          features = ncol(abundance),
          stringsAsFactors = FALSE
        )
      )
      next
    }
    matched_tax <- .match_taxonomy_to_counts(abundance, tax[[path]], allow_partial = TRUE)
    if (is.null(matched_tax) || all(is.na(.taxonomy_lineages(matched_tax)))) next
    out[[path]] <- list(
      counts = abundance,
      taxonomy = matched_tax,
      input_type = input_type
    )
  }
  attr(out, "excluded") <- excluded
  out
}

.run_study_functional <- function(parsed, study_dir, mode, tools = NULL) {
  if (mode == "off") return(.empty_functional_result(FALSE))
  .restore_functional_archive(study_dir)
  functional_dir <- file.path(study_dir, "functional")
  dir.create(functional_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(tools)) tools <- .discover_functional_tools()
  out <- .empty_functional_result(TRUE)
  manifest <- .functional_manifest_template()

  picrust_sources <- .discover_picrust_sources(study_dir)
  eligible_picrust <- 0L
  if (!nzchar(tools$picrust2)) {
    manifest <- rbind(manifest, .manifest_row(
      "picrust2", "study", "skipped", "picrust2_pipeline.py was not found on PATH.", study_dir
    ))
  } else for (source in picrust_sources) {
    if (!identical(source$modality, "amplicon")) {
      manifest <- rbind(manifest, .manifest_row(
        "picrust2", source$id, "skipped",
        "PICRUSt2 is currently supported only for amplicon ASV inputs.",
        source$counts_file
      ))
      next
    }
    loaded <- tryCatch(.load_picrust_source(source), error = function(e) e)
    branch <- source$id
    source_label <- source$counts_file
    if (inherits(loaded, "error")) {
      manifest <- rbind(manifest, .manifest_row(
        "picrust2", branch, "failed", conditionMessage(loaded), source_label
      ))
      next
    }
    if (!isTRUE(loaded$prokaryotic)) {
      manifest <- rbind(manifest, .manifest_row(
        "picrust2", branch, "skipped", "Taxonomy did not identify this branch as bacterial/archaeal 16S.",
        source_label, nrow(loaded$counts), ncol(loaded$counts)
      ))
      next
    }
    eligible_picrust <- eligible_picrust + 1L
    fit <- tryCatch(
      .run_picrust2(loaded$counts, loaded$sequences, branch, source_label,
                    functional_dir, mode, tools, loaded$input_hash,
                    taxonomy = loaded$taxonomy),
      error = function(e) list(result = NULL, manifest = .manifest_row(
        "picrust2", branch, "failed", conditionMessage(e), source_label,
        nrow(loaded$counts), ncol(loaded$counts)
      ))
    )
    manifest <- rbind(manifest, fit$manifest)
    if (!is.null(fit$result)) out$picrust2 <- .set_nested(out$picrust2, branch, fit$result)
    rm(loaded)
    gc(FALSE)
  }
  if (nzchar(tools$picrust2) && !eligible_picrust && !length(picrust_sources)) {
    manifest <- rbind(manifest, .manifest_row(
      "picrust2", "study", "skipped", "No ASV-level count and sequence input was found.", study_dir
    ))
  }

  fap_ready <- nzchar(tools$faprotax) && nzchar(tools$faprotax_db)
  fap_pairs <- if (fap_ready) .functional_faprotax_pairs(parsed) else list()
  fap_excluded <- attr(fap_pairs, "excluded", exact = TRUE)
  if (is.null(fap_excluded)) fap_excluded <- .functional_faprotax_exclusion_template()
  if (!fap_ready) {
    reason <- if (!nzchar(tools$faprotax)) {
      "collapse_table.py was not found on PATH."
    } else {
      "FAPROTAX.txt was not found beside the installed FAPROTAX executable."
    }
    manifest <- rbind(manifest, .manifest_row("faprotax", "study", "skipped", reason, study_dir))
  } else {
    if (nrow(fap_excluded)) {
      for (index in seq_len(nrow(fap_excluded))) {
        manifest <- rbind(manifest, .manifest_row(
          "faprotax", fap_excluded$branch[[index]], "skipped",
          fap_excluded$reason[[index]], fap_excluded$source[[index]],
          fap_excluded$samples[[index]], fap_excluded$features[[index]]
        ))
      }
    }
    for (branch in names(fap_pairs)) {
    pair <- fap_pairs[[branch]]
    fit <- tryCatch(
      .run_faprotax(pair$counts, pair$taxonomy, branch,
                    paste0(pair$input_type, ":", branch),
                    functional_dir, mode, tools),
      error = function(e) list(result = NULL, manifest = .manifest_row(
        "faprotax", branch, "failed", conditionMessage(e), branch,
        nrow(pair$counts), ncol(pair$counts)
      ))
    )
    manifest <- rbind(manifest, fit$manifest)
    if (!is.null(fit$result)) {
      out$faprotax <- .set_nested(out$faprotax, strsplit(branch, "/", fixed = TRUE)[[1L]], fit$result)
    }
    }
  }
  if (fap_ready && !length(fap_pairs) && !nrow(fap_excluded)) {
    manifest <- rbind(manifest, .manifest_row(
      "faprotax", "study", "skipped", "No matched bacterial/archaeal count and taxonomy branch was found.",
      study_dir
    ))
  }

  out$manifest <- manifest
  jsonlite::write_json(manifest, file.path(functional_dir, "manifest.json"),
                       dataframe = "rows", auto_unbox = TRUE, pretty = TRUE, na = "null")
  out
}
