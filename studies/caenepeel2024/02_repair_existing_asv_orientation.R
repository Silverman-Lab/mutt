suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
study_dir <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
output_dir <- file.path(study_dir, "orientation_repair")
staging_base <- paste0(output_dir, ".staging")

if (dir.exists(output_dir)) {
  stop("Repair output directory already exists; refusing to overwrite it.")
}
staging_dir <- staging_base
staging_index <- 1L
while (dir.exists(staging_dir)) {
  staging_index <- staging_index + 1L
  staging_dir <- paste0(staging_base, staging_index)
}
dir.create(staging_dir, recursive = TRUE)

threads <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(threads) || threads < 1L) threads <- 1L
threads <- min(8L, as.integer(threads))
seed <- 100L

paths <- list(
  counts_run = file.path(study_dir, "asv_counts_by_run.tsv"),
  counts_sample = file.path(study_dir, "asv_counts_by_sample_accession.tsv"),
  sequences = file.path(study_dir, "asv_sequences.tsv"),
  taxonomy_rdp16 = file.path(study_dir, "taxonomy_rdp16_toGenus.tsv"),
  taxonomy_rdp19 = file.path(study_dir, "taxonomy_rdp19_toGenus.tsv"),
  ref_rdp16 = file.path(dirname(study_dir), "helperdata", "rdp_train_set_16.fa.gz"),
  ref_rdp19 = file.path(dirname(study_dir), "helperdata", "rdp_19_toGenus_trainset.fa.gz"),
  orientation_scores = file.path(study_dir, "02_dada2_orientation_scores.cpp")
)

resolve_archived_input <- function(path) {
  if (file.exists(path)) return(path)
  archive <- paste0(path, ".zip")
  if (file.exists(archive)) return(archive)
  ""
}

zip_member <- function(path) {
  members <- unzip(path, list = TRUE)$Name
  members <- members[!grepl("/$", members)]
  expected <- sub("[.]zip$", "", basename(path))
  matching <- members[basename(members) == expected]
  if (length(matching) == 1L) return(matching)
  if (length(members) == 1L) return(members)
  stop("Archive must contain exactly one identifiable file: ", path)
}

read_archived_delim <- function(path, ...) {
  if (!grepl("[.]zip$", path)) return(read.delim(path, ...))
  connection <- unz(path, zip_member(path), open = "r")
  on.exit(if (isOpen(connection)) close(connection), add = TRUE)
  read.delim(connection, ...)
}

data_inputs <- c("counts_run", "counts_sample", "sequences", "taxonomy_rdp16", "taxonomy_rdp19")
paths[data_inputs] <- lapply(paths[data_inputs], resolve_archived_input)
missing_paths <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0L) {
  stop("Missing required input files:\n", paste(missing_paths, collapse = "\n"))
}

message("[repair] Study directory: ", study_dir)
message("[repair] Staging directory: ", staging_dir)
message("[repair] DADA2: ", as.character(packageVersion("dada2")))
message("[repair] Threads: ", threads)
message("[repair] Seed: ", seed)

read_counts <- function(path) {
  x <- read_archived_delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"SampleID" %in% names(x)) stop("Missing SampleID column: ", path)
  sample_ids <- x$SampleID
  if (any(is.na(sample_ids) | sample_ids == "") || anyDuplicated(sample_ids)) {
    stop("Sample IDs must be unique and nonblank: ", path)
  }
  x$SampleID <- NULL
  m <- as.matrix(x)
  suppressWarnings(storage.mode(m) <- "numeric")
  if (any(!is.finite(m)) || any(m < 0) || any(m != floor(m))) {
    stop("Counts must be finite nonnegative integers: ", path)
  }
  rownames(m) <- sample_ids
  m
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

reverse_complement <- function(x) {
  vapply(
    strsplit(toupper(x), "", fixed = TRUE),
    function(z) paste(chartr("ACGT", "TGCA", rev(z)), collapse = ""),
    character(1)
  )
}

prepare_reference <- function(reference) {
  refsr <- ShortRead::readFasta(reference)
  reference_sequences <- as.character(ShortRead::sread(refsr))
  keep <- nchar(reference_sequences) >= 20L
  reference_sequences <- reference_sequences[keep]
  taxonomy <- trimws(as.character(ShortRead::id(refsr))[keep])
  taxonomy_depth <- lengths(strsplit(taxonomy, ";", fixed = TRUE))
  target_depth <- max(taxonomy_depth)
  for (i in which(taxonomy_depth < target_depth)) {
    taxonomy[[i]] <- paste0(
      taxonomy[[i]],
      paste(rep("_DADA2_UNSPECIFIED;", target_depth - taxonomy_depth[[i]]), collapse = "")
    )
  }
  list(
    sequences = reference_sequences,
    reference_to_genus = match(taxonomy, unique(taxonomy))
  )
}

classify_orientation <- function(sequences, reference, label) {
  message("[repair] ", label, ": preparing DADA2 reference model")
  prepared <- prepare_reference(reference)
  message("[repair] ", label, ": comparing deterministic forward and reverse log scores")
  dada2_orientation_scores(
    sequences,
    reverse_complement(sequences),
    prepared$sequences,
    prepared$reference_to_genus
  )
}

collapse_counts <- function(counts, group_ids, ordered_ids) {
  collapsed <- t(rowsum(t(counts), group = group_ids, reorder = FALSE))
  collapsed[, ordered_ids, drop = FALSE]
}

format_taxonomy <- function(fit, asv_ids, sequences) {
  tax <- as.data.frame(fit$tax, stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(
    ASV = asv_ids,
    Sequence = sequences,
    tax,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

format_bootstraps <- function(fit, asv_ids) {
  data.frame(
    ASV = asv_ids,
    as.data.frame(fit$boot, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

message("[repair] Reading legacy ASV exports")
counts_run <- read_counts(paths$counts_run)
counts_sample <- read_counts(paths$counts_sample)
sequence_table <- read_archived_delim(
  paths$sequences, check.names = FALSE, stringsAsFactors = FALSE
)

if (!all(c("ASV", "Sequence") %in% names(sequence_table))) {
  stop("asv_sequences.tsv must contain ASV and Sequence columns.")
}
if (anyDuplicated(sequence_table$ASV) || anyDuplicated(toupper(sequence_table$Sequence))) {
  stop("Legacy sequence export contains duplicate ASV IDs or exact duplicate sequences.")
}
if (!identical(colnames(counts_run), sequence_table$ASV) ||
    !identical(colnames(counts_sample), sequence_table$ASV)) {
  stop("Count-table ASV columns do not exactly match asv_sequences.tsv.")
}

original_sequences <- toupper(sequence_table$Sequence)
if (any(!grepl("^[ACGT]+$", original_sequences))) {
  stop("Legacy ASV sequences contain characters outside A/C/G/T.")
}

message("[repair] Compiling deterministic DADA2 orientation scorer")
if (!nzchar(Sys.which("x86_64-conda-linux-gnu-c++")) && nzchar(Sys.which("g++"))) {
  makevars <- tempfile("mutt-makevars-")
  writeLines(
    c(
      "CXX = /usr/bin/g++",
      "CXX11 = /usr/bin/g++",
      "CXX14 = /usr/bin/g++",
      "CXX17 = /usr/bin/g++",
      "CXX20 = /usr/bin/g++"
    ),
    makevars
  )
  Sys.setenv(R_MAKEVARS_USER = makevars)
}
Rcpp::sourceCpp(paths$orientation_scores, rebuild = FALSE, showOutput = FALSE)

message("[repair] Determining orientation with RDP19 and RDP16")
orientation_rdp19 <- classify_orientation(original_sequences, paths$ref_rdp19, "RDP19")
orientation_rdp16 <- classify_orientation(original_sequences, paths$ref_rdp16, "RDP16")

decisive19 <- orientation_rdp19$decision %in% c("forward", "reverse")
decisive16 <- orientation_rdp16$decision %in% c("forward", "reverse")
conflict <- decisive19 & decisive16 &
  orientation_rdp19$decision != orientation_rdp16$decision

decision <- rep(NA_character_, length(original_sequences))
agree <- decisive19 & decisive16 & !conflict
decision[agree] <- orientation_rdp19$decision[agree]
decision[decisive19 & !decisive16] <- orientation_rdp19$decision[decisive19 & !decisive16]
decision[decisive16 & !decisive19] <- orientation_rdp16$decision[decisive16 & !decisive19]

palindrome <- original_sequences == reverse_complement(original_sequences)
decision[is.na(decision) & palindrome] <- "palindrome"
unresolved <- is.na(decision)

orientation_audit <- data.frame(
  original_asv = sequence_table$ASV,
  original_sequence = original_sequences,
  rdp19_decision = orientation_rdp19$decision,
  rdp16_decision = orientation_rdp16$decision,
  rdp19_forward_log_score = orientation_rdp19$forward_log_score,
  rdp19_reverse_log_score = orientation_rdp19$reverse_log_score,
  rdp19_reverse_minus_forward = orientation_rdp19$score_difference_reverse_minus_forward,
  rdp16_forward_log_score = orientation_rdp16$forward_log_score,
  rdp16_reverse_log_score = orientation_rdp16$reverse_log_score,
  rdp16_reverse_minus_forward = orientation_rdp16$score_difference_reverse_minus_forward,
  orientation_conflict = conflict,
  final_decision = ifelse(unresolved, "unresolved", decision),
  stringsAsFactors = FALSE
)
write_tsv(orientation_audit, file.path(staging_dir, "orientation_audit.tsv"))

if (any(conflict) || any(unresolved)) {
  stop(
    "Orientation audit stopped without changing counts: conflicts=", sum(conflict),
    ", unresolved=", sum(unresolved), ". Inspect orientation_audit.tsv."
  )
}

oriented_sequences <- original_sequences
reverse_rows <- decision == "reverse"
oriented_sequences[reverse_rows] <- reverse_complement(original_sequences[reverse_rows])

original_totals <- colSums(counts_run)
sequence_totals <- rowsum(
  matrix(original_totals, ncol = 1L),
  group = oriented_sequences,
  reorder = FALSE
)[, 1L]
ordered_sequences <- names(sort(sequence_totals, decreasing = TRUE))
ordered_sequences <- ordered_sequences[order(-sequence_totals[ordered_sequences], ordered_sequences)]
new_asv_ids <- paste0("ASV", seq_along(ordered_sequences))
new_id_by_sequence <- setNames(new_asv_ids, ordered_sequences)
group_ids <- unname(new_id_by_sequence[oriented_sequences])

counts_run_repaired <- collapse_counts(counts_run, group_ids, new_asv_ids)
counts_sample_repaired <- collapse_counts(counts_sample, group_ids, new_asv_ids)

if (!identical(rowSums(counts_run_repaired), rowSums(counts_run)) ||
    !identical(rowSums(counts_sample_repaired), rowSums(counts_sample))) {
  stop("Read totals were not conserved during reverse-complement collapse.")
}

mapping <- orientation_audit
mapping$oriented_sequence <- oriented_sequences
mapping$repaired_asv <- group_ids
mapping$original_total_reads <- as.numeric(original_totals)
mapping <- mapping[order(match(mapping$repaired_asv, new_asv_ids), mapping$original_asv), ]
write_tsv(mapping, file.path(staging_dir, "asv_orientation_map.tsv"))

message("[repair] Assigning RDP16 taxonomy to repaired ASVs")
set.seed(seed)
tax16 <- assignTaxonomy(
  ordered_sequences,
  paths$ref_rdp16,
  minBoot = 50,
  tryRC = FALSE,
  outputBootstraps = TRUE,
  multithread = threads
)

message("[repair] Assigning RDP19 taxonomy to repaired ASVs")
set.seed(seed)
tax19 <- assignTaxonomy(
  ordered_sequences,
  paths$ref_rdp19,
  minBoot = 50,
  tryRC = FALSE,
  outputBootstraps = TRUE,
  multithread = threads
)

repaired_reverse_match <- match(
  reverse_complement(ordered_sequences),
  ordered_sequences,
  nomatch = 0L
)
repaired_reverse_pairs <- which(repaired_reverse_match > seq_along(ordered_sequences))
if (length(repaired_reverse_pairs) > 0L) {
  stop("Repaired ASV set still contains exact reverse-complement pairs.")
}

counts_run_out <- data.frame(
  SampleID = rownames(counts_run_repaired),
  counts_run_repaired,
  check.names = FALSE
)
counts_sample_out <- data.frame(
  SampleID = rownames(counts_sample_repaired),
  counts_sample_repaired,
  check.names = FALSE
)
sequences_out <- data.frame(
  ASV = new_asv_ids,
  Sequence = ordered_sequences,
  Length = nchar(ordered_sequences),
  stringsAsFactors = FALSE
)

write_tsv(counts_run_out, file.path(staging_dir, "asv_counts_by_run.tsv"))
write_tsv(counts_sample_out, file.path(staging_dir, "asv_counts_by_sample_accession.tsv"))
write_tsv(sequences_out, file.path(staging_dir, "asv_sequences.tsv"))
write_tsv(
  format_taxonomy(tax16, new_asv_ids, ordered_sequences),
  file.path(staging_dir, "taxonomy_rdp16_toGenus.tsv")
)
write_tsv(
  format_taxonomy(tax19, new_asv_ids, ordered_sequences),
  file.path(staging_dir, "taxonomy_rdp19_toGenus.tsv")
)
write_tsv(
  format_bootstraps(tax16, new_asv_ids),
  file.path(staging_dir, "taxonomy_rdp16_bootstraps.tsv")
)
write_tsv(
  format_bootstraps(tax19, new_asv_ids),
  file.path(staging_dir, "taxonomy_rdp19_bootstraps.tsv")
)

seqtab_sequence <- counts_run_repaired
colnames(seqtab_sequence) <- ordered_sequences
saveRDS(counts_run_repaired, file.path(staging_dir, "seqtab_nochim_asv_columns.rds"))
saveRDS(seqtab_sequence, file.path(staging_dir, "seqtab_nochim_sequence_columns.rds"))
saveRDS(tax16$tax, file.path(staging_dir, "taxonomy_rdp16_toGenus_sequence_rownames.rds"))
saveRDS(tax19$tax, file.path(staging_dir, "taxonomy_rdp19_toGenus_sequence_rownames.rds"))

legacy_reverse_match <- match(
  reverse_complement(original_sequences),
  original_sequences,
  nomatch = 0L
)
legacy_reverse_pairs <- which(legacy_reverse_match > seq_along(original_sequences))

qc <- data.frame(
  asv_count = length(ordered_sequences),
  exact_reverse_complement_pair_count = length(repaired_reverse_pairs),
  all_asvs_forward_oriented = length(repaired_reverse_pairs) == 0L &&
    !any(conflict) && !any(unresolved),
  original_asv_count = length(original_sequences),
  original_reverse_complement_pair_count = length(legacy_reverse_pairs),
  original_total_reads = sum(counts_run),
  repaired_total_reads = sum(counts_run_repaired),
  sequences_reverse_complemented = sum(reverse_rows),
  rdp_orientation_conflicts = sum(conflict),
  unresolved_orientations = sum(unresolved),
  orientation_method = "DADA2 tryRC deterministic log-score comparison reproduced with RDP16 and RDP19",
  dada2_version = as.character(packageVersion("dada2")),
  seed = seed,
  stringsAsFactors = FALSE
)
write_tsv(qc, file.path(staging_dir, "asv_orientation_qc.tsv"))

manifest <- data.frame(
  file = basename(unlist(paths)[seq_len(5L)]),
  md5 = unname(tools::md5sum(unlist(paths)[seq_len(5L)])),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(staging_dir, "source_file_manifest.tsv"))

if (!file.rename(staging_dir, output_dir)) {
  stop("Repair succeeded but staging directory could not be renamed to the final output directory.")
}

message("[repair] Complete: ", output_dir)
message("[repair] ASVs: ", length(original_sequences), " -> ", length(ordered_sequences))
message("[repair] Reads conserved: ", sum(counts_run_repaired))
