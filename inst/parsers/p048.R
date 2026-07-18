parse_2024_caenepeel_gastroenterology_IBD_CD_flow <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ",
            paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }

    if (!is.logical(raw) || length(raw) != 1) {
        stop("`raw` must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(align) || length(align) != 1) {
        stop("`align` must be a single logical value (TRUE or FALSE)")
    }

    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2024_caenepeel_gastroenterology_IBD_CD_flow")

    # ----- File paths -----
    metadata_zip <- file.path(local, "SraRunTable (40).csv.zip")
    derived <- file.path(local, "orientation_repair")

    counts_by_run_tsv <- file.path(derived, "asv_counts_by_run.tsv")
    counts_by_sample_tsv <- file.path(derived, "asv_counts_by_sample_accession.tsv")
    asv_sequences_tsv <- file.path(derived, "asv_sequences.tsv")
    orientation_qc_tsv <- file.path(derived, "asv_orientation_qc.tsv")

    tax_rdp16_tsv <- file.path(derived, "taxonomy_rdp16_toGenus.tsv")
    tax_rdp19_tsv <- file.path(derived, "taxonomy_rdp19_toGenus.tsv")

    # Authors stated scale was in supplementary table 1, but did not publish it.
    scale <- NA

    counts_original <- NA
    tax_original <- NA
    proportions_original <- NA

    # -------------------------------------------------------------------------
    # Local helpers
    # -------------------------------------------------------------------------

    .check_exists <- function(paths) {
        missing <- paths[!vapply(
            paths,
            function(path) nzchar(.resolve_study_file(path, must_work = FALSE)),
            logical(1)
        )]
        if (length(missing) > 0) {
            stop(
                "Missing required file(s):\n",
                paste(missing, collapse = "\n")
            )
        }
    }

    .read_counts_tsv <- function(path) {
        x <- .read_study_delim(
            path,
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            stringsAsFactors = FALSE,
            comment.char = "",
            quote = ""
        )

        if (!"SampleID" %in% names(x)) {
            stop("Expected a `SampleID` column in: ", path)
        }

        sample_ids <- x$SampleID
        x$SampleID <- NULL

        if (any(is.na(sample_ids) | sample_ids == "") || anyDuplicated(sample_ids)) {
            stop("Count table requires unique, nonblank sample IDs: ", path)
        }

        x <- as.data.frame(
            lapply(x, as.numeric),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        count_matrix <- as.matrix(x)
        if (any(!is.finite(count_matrix))) {
            stop("Count table contains missing or non-finite values: ", path)
        }
        if (any(count_matrix < 0) || any(count_matrix != floor(count_matrix))) {
            stop("Count table must contain nonnegative integer counts: ", path)
        }

        rownames(x) <- sample_ids
        x
    }

    .read_asv_sequences <- function(path) {
        x <- .read_study_delim(
            path,
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            stringsAsFactors = FALSE,
            comment.char = "",
            quote = ""
        )

        required_cols <- c("ASV", "Sequence")
        missing_cols <- setdiff(required_cols, names(x))
        if (length(missing_cols) > 0) {
            stop(
                "`asv_sequences.tsv` is missing required column(s): ",
                paste(missing_cols, collapse = ", ")
            )
        }

        if (anyDuplicated(x$ASV)) {
            stop("`asv_sequences.tsv` contains duplicated ASV IDs.")
        }

        if (any(is.na(x$ASV) | x$ASV == "")) {
            stop("`asv_sequences.tsv` contains missing or blank ASV IDs.")
        }

        x$Sequence <- toupper(x$Sequence)
        if (any(is.na(x$Sequence) | x$Sequence == "")) {
            stop("`asv_sequences.tsv` contains missing or blank sequences.")
        }

        if (any(!grepl("^[ACGT]+$", x$Sequence))) {
            stop("`asv_sequences.tsv` contains sequences outside the A/C/G/T alphabet.")
        }
        if (anyDuplicated(x$Sequence)) {
            stop("`asv_sequences.tsv` contains duplicated sequence strings.")
        }

        reverse_complement <- function(sequence) {
            vapply(
                strsplit(sequence, "", fixed = TRUE),
                function(z) paste(chartr("ACGT", "TGCA", rev(z)), collapse = ""),
                character(1)
            )
        }
        reverse_match <- match(reverse_complement(x$Sequence), x$Sequence, nomatch = 0L)
        reverse_pairs <- which(reverse_match > seq_along(reverse_match))
        if (length(reverse_pairs) > 0) {
            stop(
                "Caenepeel ASV orientation validation failed: ",
                length(reverse_pairs),
                " exact reverse-complement ASV pair(s) remain. Rebuild the derived ",
                "counts, sequences, and taxonomy with `01_run_PRJEB71738_dada2_rdp.sbatch`."
            )
        }

        setNames(x$Sequence, x$ASV)
    }

    .read_tax_tsv_as_taxa <- function(path, seq_by_asv) {
        x <- .read_study_delim(
            path,
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            stringsAsFactors = FALSE,
            comment.char = "",
            quote = ""
        )

        required_cols <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        missing_cols <- setdiff(required_cols, names(x))
        if (length(missing_cols) > 0) {
            stop(
                "Taxonomy file is missing required column(s) in ",
                path,
                ": ",
                paste(missing_cols, collapse = ", ")
            )
        }

        if (anyDuplicated(x$ASV)) {
            stop("Taxonomy file contains duplicated ASV IDs: ", path)
        }

        if (!setequal(x$ASV, names(seq_by_asv))) {
            stop("Taxonomy and ASV sequence files do not contain the same ASV IDs: ", path)
        }

        if ("Sequence" %in% names(x)) {
            expected_sequence <- unname(seq_by_asv[x$ASV])
            observed_sequence <- toupper(x$Sequence)
            mismatch <- is.na(observed_sequence) | observed_sequence != expected_sequence
            if (any(mismatch)) {
                stop(sum(mismatch), " taxonomy sequence(s) disagree with `asv_sequences.tsv`: ", path)
            }
        }

        ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        asv_ids <- x$ASV

        tax_rank <- as.data.frame(
            x[, ranks, drop = FALSE],
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        rownames(tax_rank) <- asv_ids

        tax_labeled <- make_taxa_label(tax_rank)

        if (nrow(tax_labeled) != length(asv_ids)) {
            stop("`make_taxa_label()` changed the number of taxonomy rows for: ", path)
        }

        # Force ASV identity back after make_taxa_label().
        # Do not assume make_taxa_label() preserved rownames.
        rownames(tax_labeled) <- asv_ids
        tax_labeled$ASV <- asv_ids

        # Canonical source for sequence strings.
        tax_labeled$Sequence <- unname(seq_by_asv[asv_ids])

        # Fallback to the taxonomy TSV sequence column only if the canonical map is incomplete.
        if ("Sequence" %in% names(x)) {
            missing_seq <- is.na(tax_labeled$Sequence) | tax_labeled$Sequence == ""
            tax_labeled$Sequence[missing_seq] <- x$Sequence[missing_seq]
        }

        missing_seq <- is.na(tax_labeled$Sequence) | tax_labeled$Sequence == ""
        if (any(missing_seq)) {
            stop(
                sum(missing_seq),
                " ASVs have missing sequences after joining taxonomy against `asv_sequences.tsv`: ",
                path
            )
        }

        tax_labeled
    }

    .collapse_counts_to_taxa <- function(counts, tax) {
        if (!is.data.frame(counts) && !is.matrix(counts)) {
            return(NA)
        }

        taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]

        missing_taxa <- is.na(taxa) | taxa == "" | taxa == "NA"
        if (any(missing_taxa)) {
            taxa[missing_taxa] <- paste0("Unclassified_", colnames(counts)[missing_taxa])
        }

        colnames(counts) <- taxa
        counts <- collapse_duplicate_columns_exact(counts)

        counts <- as.data.frame(
            lapply(counts, as.numeric),
            row.names = rownames(counts),
            check.names = FALSE
        )

        counts
    }

    .proportions <- function(counts) {
        if (!is.data.frame(counts) && !is.matrix(counts)) {
            return(NA)
        }

        m <- as.matrix(counts)
        storage.mode(m) <- "numeric"

        rs <- rowSums(m, na.rm = TRUE)
        p <- sweep(m, 1, rs, "/")
        p[!is.finite(p)] <- 0

        as.data.frame(
            p,
            row.names = rownames(counts),
            check.names = FALSE
        )
    }

    .fill_if_table <- function(x) {
        if (is.data.frame(x) || is.matrix(x)) {
            return(fill_na_zero_numeric(x))
        }
        x
    }

    # -------------------------------------------------------------------------
    # Check files
    # -------------------------------------------------------------------------

    if (!nzchar(.resolve_study_file(orientation_qc_tsv, must_work = FALSE))) {
        stop(
            "The Caenepeel ASV exports predate orientation normalization. Rebuild the ",
            "derived counts, sequences, and taxonomy with either ",
            "`02_repair_existing_asv_orientation.R` or the raw-read workflow before parsing this study."
        )
    }

    .check_exists(c(
        metadata_zip,
        counts_by_run_tsv,
        counts_by_sample_tsv,
        asv_sequences_tsv,
        orientation_qc_tsv,
        tax_rdp16_tsv,
        tax_rdp19_tsv
    ))

    orientation_qc <- .read_study_delim(
        orientation_qc_tsv,
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    required_orientation_qc <- c(
        "asv_count", "exact_reverse_complement_pair_count", "all_asvs_forward_oriented"
    )
    if (nrow(orientation_qc) != 1L ||
        !all(required_orientation_qc %in% names(orientation_qc)) ||
        orientation_qc$exact_reverse_complement_pair_count[[1L]] != 0 ||
        !isTRUE(orientation_qc$all_asvs_forward_oriented[[1L]])) {
        stop("`asv_orientation_qc.tsv` does not certify a forward-oriented ASV table.")
    }

    # -------------------------------------------------------------------------
    # Metadata
    # -------------------------------------------------------------------------

    metadata <- read_zipped_table(metadata_zip, row.names = NULL)

    if (!"Run" %in% names(metadata)) {
        stop("Metadata file must contain a `Run` column.")
    }

    names(metadata)[names(metadata) == "Run"] <- "Accession"

    # rename_and_align expects `Accession` to match count-table rownames
    # and `by_col` to be the desired sample-level ID.
    if (!"Sample" %in% names(metadata)) {
        if ("BioSample" %in% names(metadata)) {
            metadata$Sample <- metadata$BioSample
        } else {
            stop("Metadata must contain either `Sample` or `BioSample`.")
        }
    }

    # -------------------------------------------------------------------------
    # Counts
    # -------------------------------------------------------------------------

    # Run-level counts are used as the main input because their rownames match metadata$Accession.
    counts_reprocessed_asv <- .read_counts_tsv(counts_by_run_tsv)

    # Read sample-accession counts to validate that the updated file exists and is parseable.
    # This object is not returned directly because rename_and_align handles run-to-sample mapping.
    counts_reprocessed_asv_by_sample <- .read_counts_tsv(counts_by_sample_tsv)

    # Basic consistency check: same ASV set in both count exports.
    if (!identical(colnames(counts_reprocessed_asv), colnames(counts_reprocessed_asv_by_sample))) {
        stop("Run-level and sample-accession-level count tables do not have identical ASV columns.")
    }

    # -------------------------------------------------------------------------
    # ASV sequence map
    # -------------------------------------------------------------------------

    seq_by_asv <- .read_asv_sequences(asv_sequences_tsv)

    if (orientation_qc$asv_count[[1L]] != length(seq_by_asv)) {
        stop("`asv_orientation_qc.tsv` ASV count does not match `asv_sequences.tsv`.")
    }

    if (!setequal(colnames(counts_reprocessed_asv), names(seq_by_asv))) {
        stop(
            "The count table and `asv_sequences.tsv` do not contain the same ASV IDs."
        )
    }

    # -------------------------------------------------------------------------
    # Taxonomy
    # -------------------------------------------------------------------------

    tax_reprocessed <- .read_tax_tsv_as_taxa(
        path = tax_rdp19_tsv,
        seq_by_asv = seq_by_asv
    )

    tax_reprocessed2 <- .read_tax_tsv_as_taxa(
        path = tax_rdp16_tsv,
        seq_by_asv = seq_by_asv
    )

    missing_rdp19_tax <- setdiff(colnames(counts_reprocessed_asv), rownames(tax_reprocessed))
    missing_rdp16_tax <- setdiff(colnames(counts_reprocessed_asv), rownames(tax_reprocessed2))

    if (length(missing_rdp19_tax) > 0) {
        stop(length(missing_rdp19_tax), " ASV column(s) are missing from the RDP19 taxonomy table.")
    }

    if (length(missing_rdp16_tax) > 0) {
        stop(length(missing_rdp16_tax), " ASV column(s) are missing from the RDP16 taxonomy table.")
    }

    # -------------------------------------------------------------------------
    # Build returned count/proportion objects
    # -------------------------------------------------------------------------

    if (raw) {
        # Raw mode returns ASV-level counts.
        counts_reprocessed <- counts_reprocessed_asv
        counts_reprocessed2 <- counts_reprocessed_asv

        proportions_reprocessed <- .proportions(counts_reprocessed)
        proportions_reprocessed2 <- .proportions(counts_reprocessed2)
    } else {
        # Convert run accessions to sample accessions.
        aligned <- rename_and_align(
            counts_reprocessed = counts_reprocessed_asv,
            metadata = metadata,
            scale = scale,
            by_col = "Sample",
            align = align,
            study_name = basename(local)
        )

        counts_reprocessed_asv <- aligned$reprocessed

        # Collapse ASVs to taxonomy labels separately for RDP19 and RDP16.
        counts_reprocessed <- .collapse_counts_to_taxa(
            counts = counts_reprocessed_asv,
            tax = tax_reprocessed
        )

        counts_reprocessed2 <- .collapse_counts_to_taxa(
            counts = counts_reprocessed_asv,
            tax = tax_reprocessed2
        )

        proportions_reprocessed <- .proportions(counts_reprocessed)
        proportions_reprocessed2 <- .proportions(counts_reprocessed2)
    }

    # -------------------------------------------------------------------------
    # Final cleanup
    # -------------------------------------------------------------------------

    if (!raw) {
        counts_original <- .fill_if_table(counts_original)
        proportions_original <- .fill_if_table(proportions_original)

        counts_reprocessed <- .fill_if_table(counts_reprocessed)
        proportions_reprocessed <- .fill_if_table(proportions_reprocessed)

        counts_reprocessed2 <- .fill_if_table(counts_reprocessed2)
        proportions_reprocessed2 <- .fill_if_table(proportions_reprocessed2)
    }

    # -------------------------------------------------------------------------
    # Return structured list
    # -------------------------------------------------------------------------

    return(list(
        counts = list(
            original = counts_original,
            reprocessed = list(
                rdp19 = counts_reprocessed,
                rdp16 = counts_reprocessed2
            )
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = list(
                rdp19 = proportions_reprocessed,
                rdp16 = proportions_reprocessed2
            )
        ),
        tax = list(
            original = tax_original,
            reprocessed = list(
                rdp19 = tax_reprocessed,
                rdp16 = tax_reprocessed2
            )
        ),
        scale = scale,
        metadata = metadata
    ))
}
