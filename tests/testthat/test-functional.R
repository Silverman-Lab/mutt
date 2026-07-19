test_that("functional mode has a small, strict public interface", {
  expect_identical(mutt:::.functional_mode(FALSE), "off")
  expect_identical(mutt:::.functional_mode(TRUE), "use")
  expect_identical(mutt:::.functional_mode("REBUILD"), "rebuild")
  expect_error(mutt:::.functional_mode("picrust2"), "FALSE, TRUE")
  expect_error(mutt:::.functional_mode(NA), "FALSE, TRUE")

  empty <- mutt:::.empty_functional_result()
  expect_named(empty, c("picrust2", "faprotax", "manifest"))
  expect_equal(nrow(empty$manifest), 0L)
})

test_that("functional process selection is internal and bounded", {
  withr::local_envvar(MUTT_FUNCTIONAL_PROCESSES = "3")
  expect_identical(mutt:::.functional_processes(), 3L)
})

test_that("functional branch modality routing is conservative", {
  expect_identical(
    mutt:::.functional_branch_modality("reprocessed/amplicon/rdp19"),
    "amplicon"
  )
  expect_identical(
    mutt:::.functional_branch_modality("reprocessed/shotgun/MetaPhlAn4"),
    "metagenomic"
  )
  expect_identical(mutt:::.functional_branch_modality("original"), "unknown")
})

test_that("tool versions can use the PICRUSt2 -v convention", {
  skip_on_os("windows")
  executable <- tempfile("mutt_version_test_")
  on.exit(unlink(executable), add = TRUE)
  writeLines(c(
    "#!/usr/bin/env bash",
    "if [[ ${1-} == -v ]]; then echo 'tool 2.6.3'; exit 0; fi",
    "exit 1"
  ), executable)
  Sys.chmod(executable, "0755")
  expect_identical(mutt:::.command_version(executable), "tool 2.6.3")
})

test_that("valid ASV sequences already stored in taxonomy are preserved", {
  tax <- data.frame(
    Kingdom = "Bacteria",
    Sequence = "ACGTRYSWKMBDHVN",
    row.names = "ASV_1"
  )
  expect_identical(mutt:::add_sequence_column(tax)$Sequence, tax$Sequence)
})

test_that("functional count validation is strict and preserves orientation", {
  counts <- matrix(
    c(1, 0, 2, 3), nrow = 2,
    dimnames = list(c("S1", "S2"), c("ACGT", "TGCA"))
  )
  expect_identical(mutt:::.normalize_count_matrix(counts, TRUE), counts)
  expect_identical(mutt:::.normalize_count_matrix(t(counts), TRUE), counts)

  bad <- counts
  bad[1, 1] <- -1
  expect_error(mutt:::.normalize_count_matrix(bad), "negative")
  bad[1, 1] <- 0.5
  expect_error(mutt:::.normalize_count_matrix(bad), "integer-valued")
})

test_that("FAPROTAX uses amplicon counts when available and proportions otherwise", {
  counts <- matrix(
    c(2L, 1L, 3L, 4L), nrow = 2L,
    dimnames = list(c("S1", "S2"), c("g_Lactobacillus", "g_Bacteroides"))
  )
  proportions <- counts / rowSums(counts)
  tax <- data.frame(
    Taxa = colnames(counts),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  with_counts <- mutt:::.functional_faprotax_pairs(list(
    counts = list(amplicon = counts),
    proportions = list(amplicon = proportions),
    tax = list(amplicon = tax)
  ))
  expect_identical(with_counts$amplicon$input_type, "counts")
  expect_equal(with_counts$amplicon$counts, counts)
  expect_identical(
    mutt:::.taxonomy_lineages(with_counts$amplicon$taxonomy),
    c("g__Lactobacillus", "g__Bacteroides")
  )

  without_counts <- mutt:::.functional_faprotax_pairs(list(
    counts = NA,
    proportions = list(amplicon = proportions),
    tax = list(amplicon = tax)
  ))
  expect_identical(without_counts$amplicon$input_type, "proportions")
  expect_equal(without_counts$amplicon$counts, proportions)
})

test_that("FAPROTAX skips metagenomic and unknown branches", {
  counts <- matrix(
    c(2L, 1L, 3L, 4L), nrow = 2L,
    dimnames = list(c("S1", "S2"), c("g_A", "g_B"))
  )
  tax <- data.frame(Taxa = colnames(counts), row.names = colnames(counts))
  pairs <- mutt:::.functional_faprotax_pairs(list(
    counts = list(
      reprocessed = list(shotgun = list(MetaPhlAn4 = counts)),
      original = counts
    ),
    proportions = NA,
    tax = list(
      reprocessed = list(shotgun = list(MetaPhlAn4 = tax)),
      original = tax
    )
  ))
  expect_length(pairs, 0L)
  excluded <- attr(pairs, "excluded", exact = TRUE)
  expect_setequal(
    excluded$branch,
    c("reprocessed/shotgun/MetaPhlAn4", "original")
  )
  expect_match(
    excluded$reason[excluded$branch == "reprocessed/shotgun/MetaPhlAn4"],
    "shotgun/metagenomic"
  )
})

test_that("PICRUSt2 rejects exact reverse-complement duplicate features", {
  expect_identical(mutt:::.reverse_complement_dna("AAGC"), "GCTT")
  expect_invisible(
    mutt:::.assert_no_reverse_complement_duplicates(c("AAGC", "ACGT"), "fixture")
  )
  expect_error(
    mutt:::.assert_no_reverse_complement_duplicates(c("AAGC", "GCTT"), "fixture"),
    "1 exact reverse-complement"
  )
})

test_that("PICRUSt2 reverse-complement repair preserves counts and taxonomy alignment", {
  counts <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    dimnames = list(c("S1", "S2"), c("AAGC", "GCTT", "CCGA"))
  )
  sequences <- setNames(colnames(counts), colnames(counts))
  taxonomy <- data.frame(
    Kingdom = rep("Bacteria", 3),
    Genus = c("Alpha", "Alpha", "Beta"),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  repaired <- mutt:::.collapse_reverse_complement_duplicates(
    counts, sequences, taxonomy
  )

  expect_equal(repaired$merged_pairs, 1L)
  expect_identical(colnames(repaired$counts), c("AAGC", "CCGA"))
  expect_equal(repaired$counts[, "AAGC"], counts[, "AAGC"] + counts[, "GCTT"])
  expect_equal(rowSums(repaired$counts), rowSums(counts))
  expect_identical(names(repaired$sequences), colnames(repaired$counts))
  expect_identical(rownames(repaired$taxonomy), colnames(repaired$counts))
  expect_identical(repaired$taxonomy$Genus, c("Alpha", "Beta"))
  expect_invisible(
    mutt:::.assert_no_reverse_complement_duplicates(repaired$sequences, "fixture")
  )
})

test_that("taxonomy is matched and formatted for FAPROTAX", {
  counts <- matrix(
    c(2, 1, 3, 4), nrow = 2,
    dimnames = list(c("S1", "S2"), c("g_Lactobacillus", "g_Bacteroides"))
  )
  tax <- data.frame(
    Kingdom = c("k__Bacteria", "Bacteria", "Bacteria"),
    Phylum = c("p__Firmicutes", "Firmicutes", "Bacteroidota"),
    Genus = c("g__Lactobacillus", "Lactobacillus", "Bacteroides"),
    Species = c("s__Lactobacillus_casei", NA, NA),
    Taxa = c("g_Lactobacillus", "g_Lactobacillus", "g_Bacteroides"),
    row.names = c("ASV1", "ASV2", "ASV3"),
    check.names = FALSE
  )

  matched <- mutt:::.match_taxonomy_to_counts(counts, tax)
  expect_true(mutt:::.is_prokaryotic_taxonomy(matched))
  lineages <- mutt:::.taxonomy_lineages(matched)
  expect_match(lineages[[1]], "^k__Bacteria;p__Firmicutes;g__Lactobacillus")
  expect_false(grepl("k__k__", lineages[[1]], fixed = TRUE))

  partial_counts <- cbind(counts, g_Unclassified = c(1, 0))
  expect_null(mutt:::.match_taxonomy_to_counts(partial_counts, tax))
  partial <- mutt:::.match_taxonomy_to_counts(
    partial_counts,
    tax,
    allow_partial = TRUE
  )
  expect_equal(rownames(partial), colnames(partial_counts))
  expect_true(all(is.na(partial["g_Unclassified", ])))
})

test_that("functional output tables retain only numeric sample columns", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "# generated output",
    "function\tS1\tS2\tdescription",
    "F1\t1\t2\tfirst",
    "F2\t3.5\t4\tsecond"
  ), path)

  observed <- mutt:::.read_function_table(path, c("S1", "S2"))
  expected <- matrix(
    c(1, 2, 3.5, 4), nrow = 2,
    dimnames = list(c("S1", "S2"), c("F1", "F2"))
  )
  expect_equal(observed, expected)
})

test_that("PICRUSt2 output restores missing input samples as diagnosed zero rows", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "function\tS1\tS3",
    "F1\t1\t3",
    "F2\t2\t4"
  ), path)

  observed <- mutt:::.read_function_table(
    path,
    c("S1", "S2", "S3"),
    allow_missing_samples = TRUE
  )
  expect_identical(rownames(observed), c("S1", "S2", "S3"))
  expect_equal(unname(observed["S2", ]), c(0, 0))
  reconciliation <- attr(observed, "sample_reconciliation", exact = TRUE)
  expect_identical(reconciliation$missing_samples, "S2")
  expect_identical(reconciliation$zero_filled_missing_samples, 1L)
})

test_that("functional output rejects unexpected sample IDs", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("function\tS1\tRENAMED", "F1\t1\t2"), path)
  expect_error(
    mutt:::.read_function_table(
      path,
      c("S1", "S2"),
      allow_missing_samples = TRUE
    ),
    "1 expected sample.*1 unexpected"
  )
})

test_that("PICRUSt2 parsing reconciles the same missing samples across outputs", {
  root <- tempfile("mutt_picrust_parse_")
  dir.create(file.path(root, "EC_metagenome_out"), recursive = TRUE)
  dir.create(file.path(root, "KO_metagenome_out"), recursive = TRUE)
  dir.create(file.path(root, "pathways_out"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  write_gzip <- function(path, lines) {
    connection <- gzfile(path, "wt")
    on.exit(close(connection), add = TRUE)
    writeLines(lines, connection)
  }
  write_gzip(
    file.path(root, "EC_metagenome_out", "pred_metagenome_unstrat.tsv.gz"),
    c("function\tS1", "EC1\t1")
  )
  write_gzip(
    file.path(root, "KO_metagenome_out", "pred_metagenome_unstrat.tsv.gz"),
    c("function\tS1", "KO1\t2")
  )
  write_gzip(
    file.path(root, "pathways_out", "path_abun_unstrat.tsv.gz"),
    c("function\tS1", "PWY1\t3")
  )
  write_gzip(
    file.path(root, "pathways_out", "path_cov_unstrat.tsv.gz"),
    c("function\tS1", "PWY1\t0.5")
  )
  contribution_header <- paste(
    c(
      "sample", "function", "taxon", "taxon_abun", "taxon_rel_abun",
      "genome_function_count", "taxon_function_abun",
      "taxon_rel_function_abun", "norm_taxon_function_contrib"
    ),
    collapse = "\t"
  )
  write_gzip(
    file.path(root, "EC_metagenome_out", "pred_metagenome_contrib.tsv.gz"),
    c(contribution_header, "S1\tEC1\tASV1\t1\t1\t1\t1\t1\t1")
  )
  write_gzip(
    file.path(root, "KO_metagenome_out", "pred_metagenome_contrib.tsv.gz"),
    c(contribution_header, "S1\tKO1\tASV1\t1\t1\t1\t1\t1\t1")
  )
  write_gzip(
    file.path(root, "pathways_out", "path_abun_contrib.tsv.gz"),
    c(contribution_header, "S1\tPWY1\tASV1\t1\t1\t1\t1\t1\t1")
  )

  counts <- matrix(
    c(2L, 0L, 1L, 3L),
    nrow = 2L,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )
  parsed <- mutt:::.parse_picrust_output(root, counts)
  expect_equal(parsed$ec["S2", "EC1"], 0)
  expect_equal(parsed$ko["S2", "KO1"], 0)
  expect_equal(parsed$metacyc_abundance["S2", "PWY1"], 0)
  expect_equal(parsed$metacyc_coverage["S2", "PWY1"], 0)
  expect_identical(parsed$qc$picrust_output_samples, 1L)
  expect_identical(parsed$qc$zero_filled_missing_samples, 1L)
  expect_identical(
    parsed$sample_reconciliation$status,
    c("retained", "zero_filled_missing")
  )
})

test_that("PICRUSt2 discovery ignores archives disabled in a parser", {
  root <- tempfile("mutt_discovery_test_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  active <- file.path(root, "ACTIVE_dada2_counts.rds.zip")
  inactive <- file.path(root, "INACTIVE_dada2_counts.rds.zip")
  file.create(active, inactive)
  writeLines(c(
    "counts <- 'ACTIVE_dada2_counts.rds.zip'",
    "# counts <- 'INACTIVE_dada2_counts.rds.zip'"
  ), file.path(root, "parse.R"))

  sources <- mutt:::.discover_picrust_sources(root)
  expect_identical(vapply(sources, `[[`, character(1), "id"), "ACTIVE")
  expect_identical(vapply(sources, `[[`, character(1), "modality"), "amplicon")
})

test_that("PICRUSt2 post-processing failures retain staging output", {
  root <- tempfile("mutt_picrust_failed_staging_")
  staging <- file.path(root, "staging")
  target <- file.path(root, "picrust2", "branch")
  dir.create(file.path(staging, "raw"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  writeLines("diagnostic", file.path(staging, "raw", "output.txt"))

  retained <- mutt:::.retain_failed_picrust_staging(staging, target)
  expect_false(dir.exists(staging))
  expect_true(file.exists(file.path(retained, "raw", "output.txt")))
  expect_match(retained, "branch[.]failed")
})

test_that("explicit ASV discovery prefers validated orientation repair outputs", {
  root <- tempfile("mutt_orientation_repair_test_")
  repaired <- file.path(root, "orientation_repair")
  dir.create(repaired, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  required <- c(
    "asv_counts_by_run.tsv",
    "asv_sequences.tsv",
    "asv_orientation_qc.tsv",
    "taxonomy_rdp19_toGenus.tsv"
  )
  file.create(file.path(root, required))
  file.create(file.path(repaired, required))

  sources <- mutt:::.discover_picrust_sources(root)
  explicit <- sources[vapply(sources, `[[`, character(1), "type") == "explicit"]
  expect_length(explicit, 1L)
  expect_identical(
    normalizePath(dirname(explicit[[1L]]$counts_file), winslash = "/"),
    normalizePath(repaired, winslash = "/")
  )
})

test_that("explicit ASV discovery and loading support one-file ZIP archives", {
  skip_if_not(nzchar(Sys.which("zip")))
  withr::local_envvar(R_ZIPCMD = Sys.which("zip"))
  root <- tempfile("mutt_zipped_asv_test_")
  repaired <- file.path(root, "orientation_repair")
  dir.create(repaired, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  write.table(
    data.frame(SampleID = "S1", ASV_1 = 4L, check.names = FALSE),
    file.path(repaired, "asv_counts_by_run.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(ASV = "ASV_1", Sequence = "AAGC"),
    file.path(repaired, "asv_sequences.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(
      asv_count = 1L,
      exact_reverse_complement_pair_count = 0L,
      all_asvs_forward_oriented = TRUE
    ),
    file.path(repaired, "asv_orientation_qc.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(ASV = "ASV_1", Kingdom = "Bacteria", Genus = "Bacteroides"),
    file.path(repaired, "taxonomy_rdp19_toGenus.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  rds_value <- matrix(1:4, nrow = 2L)
  saveRDS(rds_value, file.path(repaired, "fixture.rds"))
  files <- list.files(repaired, full.names = TRUE)
  for (file in files) {
    withr::with_dir(
      repaired,
      utils::zip(paste0(basename(file), ".zip"), basename(file), flags = "-q -X")
    )
    unlink(file)
  }

  expect_identical(
    mutt:::.read_study_rds(file.path(repaired, "fixture.rds")),
    rds_value
  )

  sources <- mutt:::.discover_picrust_sources(root)
  explicit <- sources[vapply(sources, `[[`, character(1), "type") == "explicit"]
  expect_length(explicit, 1L)
  expect_true(all(grepl("[.]zip$", unlist(explicit[[1L]][c(
    "counts_file", "sequences_file", "tax_file", "orientation_qc_file"
  )]))))
  loaded <- mutt:::.load_picrust_source(explicit[[1L]])
  expect_identical(unname(loaded$counts), matrix(4, nrow = 1L))
  expect_identical(unname(loaded$sequences), "AAGC")
})

test_that("functional cache archives restore atomically from one safe root", {
  skip_if_not(nzchar(Sys.which("zip")))
  withr::local_envvar(R_ZIPCMD = Sys.which("zip"))
  study <- tempfile("mutt_functional_archive_")
  dir.create(file.path(study, "functional", "picrust2"), recursive = TRUE)
  on.exit(unlink(study, recursive = TRUE), add = TRUE)
  writeLines("cache", file.path(study, "functional", "picrust2", "marker.txt"))
  withr::with_dir(
    study,
    utils::zip("functional.zip", "functional", flags = "-q -X -r")
  )
  unlink(file.path(study, "functional"), recursive = TRUE)

  expect_true(mutt:::.restore_functional_archive(study))
  expect_identical(
    readLines(file.path(study, "functional", "picrust2", "marker.txt")),
    "cache"
  )
  expect_false(mutt:::.restore_functional_archive(study))
})

test_that("functional cache archives reject files outside functional root", {
  skip_if_not(nzchar(Sys.which("zip")))
  withr::local_envvar(R_ZIPCMD = Sys.which("zip"))
  study <- tempfile("mutt_unsafe_functional_archive_")
  dir.create(study)
  on.exit(unlink(study, recursive = TRUE), add = TRUE)
  writeLines("unexpected", file.path(study, "unexpected.txt"))
  withr::with_dir(
    study,
    utils::zip("functional.zip", "unexpected.txt", flags = "-q -X")
  )
  expect_error(mutt:::.restore_functional_archive(study), "unsafe or unexpected")
  expect_false(dir.exists(file.path(study, "functional")))
})

test_that("PICRUSt2 cache is reusable and a failed rebuild preserves it", {
  skip_on_os("windows")
  root <- tempfile("mutt_picrust_test_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  executable <- file.path(root, "fake-picrust2")
  biom_executable <- file.path(root, "biom")
  hmmalign_executable <- file.path(root, "hmmalign")
  reference_msa <- file.path(root, "reference.fna")
  reference_hmm <- file.path(root, "reference.hmm")
  file.create(reference_msa, reference_hmm)
  writeLines(c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "if [[ ${1-} == --version ]]; then echo 'biom 2.1.17'; exit 0; fi",
    "if [[ ${1-} == validate-table ]]; then exit 0; fi",
    "if [[ ${1-} == convert ]]; then",
    "  in=''",
    "  out=''",
    "  while [[ $# -gt 0 ]]; do",
    "    if [[ $1 == -i ]]; then in=$2; shift 2",
    "    elif [[ $1 == -o ]]; then out=$2; shift 2",
    "    else shift; fi",
    "  done",
    "  cp \"$in\" \"$out\"",
    "  exit 0",
    "fi",
    "exit 1"
  ), biom_executable)
  Sys.chmod(biom_executable, "0755")
  writeLines(c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "out=''",
    "input=${@: -1}",
    "while [[ $# -gt 0 ]]; do",
    "  if [[ $1 == -o ]]; then out=$2; shift 2; else shift; fi",
    "done",
    "{",
    "  echo '# STOCKHOLM 1.0'",
    "  awk '/^>/{id=substr($0,2); next} {seq=$0; if (id ~ /^forward_/) seq=substr(seq,1,2); print id \"\\t\" seq}' \"$input\"",
    "  echo '//'",
    "} > \"$out\""
  ), hmmalign_executable)
  Sys.chmod(hmmalign_executable, "0755")
  writeLines(c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "if [[ ${1-} == --version ]]; then echo 'fake-picrust2 1'; exit 0; fi",
    "printf '%s\\n' \"$*\" > \"$(dirname \"$0\")/picrust_args.txt\"",
    "out=''",
    "while [[ $# -gt 0 ]]; do if [[ $1 == -o ]]; then out=$2; shift 2; else shift; fi; done",
    "mkdir -p \"$out/EC_metagenome_out\" \"$out/KO_metagenome_out\" \"$out/pathways_out\"",
    "printf 'function\\tS1\\tS2\\nEC1\\t1\\t2\\n' | gzip -c > \"$out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz\"",
    "printf 'sample\\tfunction\\ttaxon\\ttaxon_abun\\ttaxon_rel_abun\\tgenome_function_count\\ttaxon_function_abun\\ttaxon_rel_function_abun\\tnorm_taxon_function_contrib\\nS1\\tEC1\\tASV1\\t1\\t0.5\\t2\\t2\\t1\\t1\\n' | gzip -c > \"$out/EC_metagenome_out/pred_metagenome_contrib.tsv.gz\"",
    "printf 'function\\tS1\\tS2\\nKO1\\t3\\t4\\n' | gzip -c > \"$out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz\"",
    "printf 'sample\\tfunction\\ttaxon\\ttaxon_abun\\ttaxon_rel_abun\\tgenome_function_count\\ttaxon_function_abun\\ttaxon_rel_function_abun\\tnorm_taxon_function_contrib\\nS1\\tKO1\\tASV1\\t1\\t0.5\\t2\\t2\\t1\\t1\\n' | gzip -c > \"$out/KO_metagenome_out/pred_metagenome_contrib.tsv.gz\"",
    "printf 'function\\tS1\\tS2\\nPWY1\\t5\\t6\\n' | gzip -c > \"$out/pathways_out/path_abun_unstrat.tsv.gz\"",
    "printf 'sample\\tfunction\\ttaxon\\ttaxon_abun\\ttaxon_rel_abun\\tgenome_function_count\\ttaxon_function_abun\\ttaxon_rel_function_abun\\tnorm_taxon_function_contrib\\nS1\\tPWY1\\tASV1\\t1\\t0.5\\t2\\t2\\t1\\t1\\n' | gzip -c > \"$out/pathways_out/path_abun_contrib.tsv.gz\"",
    "printf 'function\\tS1\\tS2\\nPWY1\\t0.5\\t0.6\\n' | gzip -c > \"$out/pathways_out/path_cov_unstrat.tsv.gz\"",
    "printf 'sequence\\tNSTI\\nASV1\\t0.1\\n' | gzip -c > \"$out/marker_predicted_and_nsti.tsv.gz\""
  ), executable)
  Sys.chmod(executable, "0755")

  counts <- matrix(
    c(1, 0, 2, 3), nrow = 2,
    dimnames = list(c("S1", "S2"), c("ASV_A", "ASV_B"))
  )
  sequences <- c(ASV_A = "AAGC", ASV_B = "CCGA")
  taxonomy <- data.frame(
    Kingdom = c("Bacteria", "Bacteria"),
    Genus = c("Alpha", "Beta"),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  dependencies <- c(
    place_seqs = executable,
    hsp = executable,
    metagenome = executable,
    pathway = executable,
    biom = biom_executable,
    hmmalign = hmmalign_executable
  )
  tools <- list(
    picrust2 = executable,
    picrust2_version = "fake-picrust2 1",
    picrust2_dependencies = dependencies,
    picrust2_references = c(
      bacteria_msa = reference_msa,
      bacteria_hmm = reference_hmm,
      archaea_msa = reference_msa,
      archaea_hmm = reference_hmm
    ),
    biom_version = "biom 2.1.17"
  )

  first <- mutt:::.run_picrust2(
    counts, sequences, "asv", "fixture", root, "use", tools,
    taxonomy = taxonomy
  )
  expect_identical(first$manifest$status, "generated")
  expect_equal(first$result$ko["S2", "KO1"], 4)
  cache_dir <- file.path(root, "picrust2", "asv")
  expect_true(file.exists(file.path(cache_dir, "input", "study_sequences.fna")))
  expect_true(file.exists(file.path(cache_dir, "input", "asv_counts.biom")))
  expect_true(file.exists(file.path(cache_dir, "input", "asv_id_map.tsv")))
  expect_false(dir.exists(file.path(cache_dir, "input", "orientation")))
  expect_false(dir.exists(file.path(cache_dir, "raw", "intermediate")))
  expect_identical(first$result$asv_mapping$original_feature_id, colnames(counts))
  expect_identical(first$result$asv_mapping$taxonomy_row_id, colnames(counts))
  expect_identical(first$result$taxonomy$original_feature_id, colnames(counts))
  expect_true(all(first$result$asv_mapping$orientation == "reverse_complement"))
  expect_equal(first$result$qc$reverse_complemented_asvs, 2L)
  expect_match(readLines(file.path(root, "picrust_args.txt")), "--coverage")
  expect_match(readLines(file.path(root, "picrust_args.txt")), "--stratified")
  expect_identical(
    names(first$result$stratified),
    c("ec", "ko", "metacyc_abundance")
  )
  contributions <- as.data.frame(first$result, type = "ec")
  expect_identical(contributions$picrust_id, "ASV1")
  expect_identical(contributions$original_feature_id, "ASV_A")
  expect_identical(contributions$Genus, "Alpha")
  expect_equal(contributions$norm_taxon_function_contrib, 1)

  cached <- mutt:::.run_picrust2(
    counts, sequences, "asv", "fixture", root, "use",
    list(picrust2 = "", picrust2_version = "fake-picrust2 1"),
    taxonomy = taxonomy
  )
  expect_identical(cached$manifest$status, "cached")
  cached_contributions <- as.data.frame(cached$result, type = "ec")
  expect_identical(cached_contributions$original_feature_id, "ASV_A")
  expect_identical(cached_contributions$Genus, "Alpha")
  expect_error(
    mutt:::.run_picrust2(
      counts, sequences, "asv", "fixture", root, "rebuild",
      list(
        picrust2 = "/bin/false",
        picrust2_version = "broken",
        picrust2_dependencies = dependencies,
        picrust2_references = c(
          bacteria_msa = reference_msa,
          bacteria_hmm = reference_hmm,
          archaea_msa = reference_msa,
          archaea_hmm = reference_hmm
        ),
        biom_version = "biom 2.1.17"
      ),
      taxonomy = taxonomy
    ),
    "PICRUSt2 failed"
  )
  expect_true(file.exists(file.path(cache_dir, "result.rds")))
  expect_true(dir.exists(paste0(cache_dir, ".failed")))

  relocated <- tempfile("mutt_picrust_relocated_")
  dir.create(file.path(relocated, "picrust2"), recursive = TRUE)
  on.exit(unlink(relocated, recursive = TRUE), add = TRUE)
  expect_true(file.copy(cache_dir, file.path(relocated, "picrust2"), recursive = TRUE))
  relocated_fit <- mutt:::.run_picrust2(
    counts, sequences, "asv", "fixture", relocated, "use",
    list(picrust2 = "", picrust2_version = "fake-picrust2 1"),
    taxonomy = taxonomy
  )
  expect_identical(relocated_fit$manifest$status, "cached")
  expect_identical(
    relocated_fit$result$provenance$output_directory,
    file.path(relocated, "picrust2", "asv")
  )
  expect_identical(
    as.data.frame(relocated_fit$result, type = "ec")$original_feature_id,
    "ASV_A"
  )
})

test_that("FAPROTAX receives taxonomy and returns sample-by-function output", {
  skip_on_os("windows")
  root <- tempfile("mutt_faprotax_test_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  executable <- file.path(root, "fake-faprotax")
  database <- file.path(root, "FAPROTAX.txt")
  writeLines("fake database", database)
  writeLines(c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "if [[ ${1-} == --version ]]; then echo 'fake-faprotax 1'; exit 0; fi",
    "out=''",
    "report=''",
    "while [[ $# -gt 0 ]]; do",
    "  if [[ $1 == -o ]]; then out=$2; shift 2",
    "  elif [[ $1 == -r ]]; then report=$2; shift 2",
    "  else shift; fi",
    "done",
    "printf '#Group\\tS1\\tS2\\nfermentation\\t3\\t4\\n' > \"$out\"",
    "printf 'fermentation: g__Lactobacillus\\n' > \"$report\""
  ), executable)
  Sys.chmod(executable, "0755")

  counts <- matrix(
    c(1, 2, 3, 4), nrow = 2,
    dimnames = list(c("S1", "S2"), c("g_Lactobacillus", "g_Bacteroides"))
  )
  taxonomy <- data.frame(
    Kingdom = c("Bacteria", "Bacteria"),
    Genus = c("Lactobacillus", "Bacteroides"),
    row.names = colnames(counts)
  )
  tools <- list(
    faprotax = executable,
    faprotax_version = "fake-faprotax 1",
    faprotax_db = database,
    faprotax_db_version = unname(tools::md5sum(database))
  )

  fit <- mutt:::.run_faprotax(
    counts, taxonomy, "reprocessed/rdp19", "fixture", root, "use", tools
  )
  expect_identical(fit$manifest$status, "generated")
  expect_equal(fit$result$abundance["S2", "fermentation"], 4)
  expect_true(nrow(fit$result$assignments) > 0)
})
