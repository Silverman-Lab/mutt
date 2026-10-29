options(width = 160)

files <- c(
  "asv_counts_by_run.tsv",
  "asv_counts_by_sample_accession.tsv",
  "asv_sequences.tsv",
  "seqtab_nochim_asv_columns.rds",
  "seqtab_nochim_sequence_columns.rds",
  "SraRunTable (40).csv.zip",
  "taxonomy_rdp16_toGenus_sequence_rownames.rds",
  "taxonomy_rdp16_toGenus.tsv",
  "taxonomy_rdp19_toGenus_sequence_rownames.rds",
  "taxonomy_rdp19_toGenus.tsv"
)

resolve_archived_file <- function(path) {
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

open_archived_file <- function(path, mode) {
  if (!grepl("[.]zip$", path)) return(file(path, open = mode))
  unz(path, zip_member(path), open = mode)
}

cat("\n==============================\n")
cat("FILES PRESENT\n")
cat("==============================\n")
resolved_files <- setNames(vapply(files, resolve_archived_file, character(1)), files)
print(file.info(resolved_files[nzchar(resolved_files)])[, c("size", "mtime")])

read_tsv <- function(path) {
  connection <- open_archived_file(path, "r")
  on.exit(if (isOpen(connection)) close(connection), add = TRUE)
  read.delim(
    connection,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
}

inspect_df <- function(path) {
  cat("\n\n==============================\n")
  cat(path, "\n")
  cat("==============================\n")

  x <- read_tsv(path)

  cat("class:", paste(class(x), collapse = ", "), "\n")
  cat("dim:", paste(dim(x), collapse = " x "), "\n")
  cat("first 12 column names:\n")
  print(head(names(x), 12))

  cat("\ncolumn classes:\n")
  print(table(vapply(x, function(z) class(z)[1], character(1))))

  cat("\nfirst 6 rows x first 8 columns:\n")
  print(utils::head(x[, seq_len(min(8, ncol(x))), drop = FALSE], 6))

  numeric_cols <- vapply(x, is.numeric, logical(1))
  if (sum(numeric_cols) > 0) {
    cat("\nnumeric column count:", sum(numeric_cols), "of", ncol(x), "\n")
    cat("numeric matrix summary, using numeric columns only:\n")
    num <- as.matrix(x[, numeric_cols, drop = FALSE])
    cat("min:", min(num, na.rm = TRUE), "\n")
    cat("median:", median(num, na.rm = TRUE), "\n")
    cat("max:", max(num, na.rm = TRUE), "\n")
    cat("zero fraction:", mean(num == 0, na.rm = TRUE), "\n")

    if (nrow(num) > 1 && ncol(num) > 1) {
      cat("\nrow sums summary:\n")
      print(summary(rowSums(num, na.rm = TRUE)))
      cat("\ncolumn sums summary:\n")
      print(summary(colSums(num, na.rm = TRUE)))
    }
  }

  invisible(x)
}

inspect_rds <- function(path) {
  cat("\n\n==============================\n")
  cat(path, "\n")
  cat("==============================\n")

  if (grepl("[.]zip$", path)) {
    temporary <- tempfile("caenepeel-rds-")
    dir.create(temporary)
    on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
    extracted <- unzip(
      path,
      files = zip_member(path),
      exdir = temporary,
      junkpaths = TRUE
    )
    if (length(extracted) != 1L || !file.exists(extracted)) {
      stop("Failed to extract RDS archive: ", path)
    }
    x <- readRDS(extracted)
  } else {
    x <- readRDS(path)
  }

  cat("class:", paste(class(x), collapse = ", "), "\n")
  cat("typeof:", typeof(x), "\n")

  if (!is.null(dim(x))) {
    cat("dim:", paste(dim(x), collapse = " x "), "\n")
  } else {
    cat("length:", length(x), "\n")
  }

  cat("\nstructure, max.level = 2:\n")
  str(x, max.level = 2)

  if (is.matrix(x) || is.data.frame(x)) {
    cat("\nfirst row names:\n")
    print(head(rownames(x), 6))

    cat("\nfirst column names:\n")
    print(head(colnames(x), 6))

    if (is.matrix(x) && is.numeric(x)) {
      cat("\nnumeric matrix summary:\n")
      cat("min:", min(x, na.rm = TRUE), "\n")
      cat("median:", median(x, na.rm = TRUE), "\n")
      cat("max:", max(x, na.rm = TRUE), "\n")
      cat("zero fraction:", mean(x == 0, na.rm = TRUE), "\n")

      cat("\nrow sums summary:\n")
      print(summary(rowSums(x, na.rm = TRUE)))

      cat("\ncolumn sums summary:\n")
      print(summary(colSums(x, na.rm = TRUE)))
    }
  }

  invisible(x)
}

cat("\n\n==============================\n")
cat("TSV INSPECTION\n")
cat("==============================\n")

tsv_objects <- list()
for (f in files[grepl("\\.tsv$", files)]) {
  resolved <- resolved_files[[f]]
  if (nzchar(resolved)) {
    tsv_objects[[f]] <- inspect_df(resolved)
  } else {
    cat("\nMissing:", f, "\n")
  }
}

cat("\n\n==============================\n")
cat("RDS INSPECTION\n")
cat("==============================\n")

rds_objects <- list()
for (f in files[grepl("\\.rds$", files)]) {
  resolved <- resolved_files[[f]]
  if (nzchar(resolved)) {
    rds_objects[[f]] <- inspect_rds(resolved)
  } else {
    cat("\nMissing:", f, "\n")
  }
}

cat("\n\n==============================\n")
cat("SRA RUN TABLE INSPECTION\n")
cat("==============================\n")

sra_zip <- "SraRunTable (40).csv.zip"
if (file.exists(sra_zip)) {
  cat("zip contents:\n")
  print(unzip(sra_zip, list = TRUE))

  sra <- read.csv(unz(sra_zip, unzip(sra_zip, list = TRUE)$Name[1]), stringsAsFactors = FALSE, check.names = FALSE)
  cat("\nSRA table dim:", paste(dim(sra), collapse = " x "), "\n")
  cat("first 20 column names:\n")
  print(head(names(sra), 20))
  cat("\nfirst 6 rows x first 10 columns:\n")
  print(head(sra[, seq_len(min(10, ncol(sra))), drop = FALSE], 6))
}

cat("\n\n==============================\n")
cat("CROSS-CHECKS\n")
cat("==============================\n")

if (all(c("seqtab_nochim_asv_columns.rds", "seqtab_nochim_sequence_columns.rds") %in% names(rds_objects))) {
  a <- rds_objects[["seqtab_nochim_asv_columns.rds"]]
  s <- rds_objects[["seqtab_nochim_sequence_columns.rds"]]

  cat("\nseqtab_nochim_asv_columns.rds vs seqtab_nochim_sequence_columns.rds\n")
  cat("same dim:", identical(dim(a), dim(s)), "\n")
  cat("same rownames:", identical(rownames(a), rownames(s)), "\n")
  cat("same numeric counts ignoring column names:", isTRUE(all.equal(unname(a), unname(s))), "\n")

  cat("\nInterpretation: if TRUE/TRUE/TRUE above, these are probably the same count table, but one has ASV IDs as columns and the other has raw DNA sequences as columns.\n")
}

if (all(c("taxonomy_rdp16_toGenus_sequence_rownames.rds", "taxonomy_rdp19_toGenus_sequence_rownames.rds") %in% names(rds_objects))) {
  tax16 <- rds_objects[["taxonomy_rdp16_toGenus_sequence_rownames.rds"]]
  tax19 <- rds_objects[["taxonomy_rdp19_toGenus_sequence_rownames.rds"]]

  cat("\nRDP16 vs RDP19 taxonomy RDS\n")
  cat("same dim:", identical(dim(tax16), dim(tax19)), "\n")
  cat("same rownames:", identical(rownames(tax16), rownames(tax19)), "\n")
  cat("same column names:", identical(colnames(tax16), colnames(tax19)), "\n")
  cat("number of exact full-row taxonomy matches:", sum(apply(tax16 == tax19, 1, all), na.rm = TRUE), "\n")
  cat("number of taxa/ASVs compared:", nrow(tax16), "\n")
}

cat("\n\nDone.\n")
