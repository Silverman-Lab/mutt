parse_2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    if (!is.logical(raw)) {
        stop("raw must be a logical value")
    }
    if (!is.logical(align)) {
        stop("align must be a logical value")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct")

    # ----- File paths -----
    repro_counts_zips <- c(
    file.path(local, "PRJNA545312_dada2_counts.rds.zip"),
    file.path(local, "PRJNA548153_dada2_counts.rds.zip"),
    file.path(local, "PRJNA606262_dada2_counts.rds.zip"),
    file.path(local, "PRJNA607574_dada2_counts.rds.zip")
    )

    repro_tax_zips <- c(
    file.path(local, "PRJNA545312_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA548153_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA606262_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA607574_dada2_taxa.rds.zip")
    )

    sra_zip      <- c(
        file.path(local, "SraRunTable (53).csv.zip"),
        file.path(local, "SraRunTable (54).csv.zip"),
        file.path(local, "SraRunTable (55).csv.zip"),
        file.path(local, "SraRunTable (56).csv.zip"),
        file.path(local, "SraRunTable (57).csv.zip")
    )

    scale_zip    <- file.path(local, "Liao2021_scale.csv.zip")
    metadata_zip <- file.path(local, "Liao_2021_metadata.csv.zip")
    counts_zip   <- file.path(local, "Liao_2021_16S.csv.zip")


    dfs <- lapply(sra_zip, function(zf) {
    df <- read_zipped_table(zf, row.names = NULL)
    df <- as.data.frame(
        lapply(df, as.character),
        stringsAsFactors = FALSE
    )
    df
    })

    sracombined <- bind_rows(dfs)
    names(sracombined)[ names(sracombined) == "Run" ] <- "Accession"

    # ---- scale and metadata -----
    scale <- read_zipped_table(scale_zip, row.names = NULL) %>%
    # 1) rename your original qPCR16S → log10_qPCR
    rename(log10_qPCR = qPCR16S) %>%
    # 2) compute the base-2 version
    mutate(log2_qPCR = log10_qPCR * log2(10)) %>%
    # 3) drop any rows where either measure is 0, ±Inf or NA
    filter(
        if_all(
        c(log10_qPCR, log2_qPCR),
        ~ !is.na(.) & is.finite(.) & . != 0
        )
    )
    metadata  <- read_zipped_table(metadata_zip, row.names = NULL)

    meta_and_scale <- metadata %>%
        full_join(scale, by = "SampleID")

    metadata2 <- sracombined %>%
        left_join(meta_and_scale, by = "Accession")

    has_qPCR <- metadata2 %>% filter(!is.na(log10_qPCR))

    # ------ original counts ------
    counts_original <- read_zipped_table(counts_zip)

    if (!is.na(counts_original)[1]) {
        original_taxa <- colnames(counts_original)

        # Create taxa mapping data frame
        tax_original <- data.frame(
        Taxa = original_taxa,
        stringsAsFactors = FALSE 
        )

        if (!raw) {
            aligned <- rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
            counts_original <- aligned$counts_original
            original_names <- colnames(counts_original)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
        }
        # proportions
        proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")
    } else {
        proportions_original <- NA
        tax_original <- NA
    }

    # Process multiple zipped RDS files
    counts_reprocessed_list <- list()
    proportions_reprocessed_list <- list()
    tax_reprocessed_list <- list()
    tax_reprocessed2_list <- list()
    counts_reprocessed2_list <- list()
    proportions_reprocessed2_list <- list() 

    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    tax_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    counts_reprocessed2 <- NA

    if (all(file.exists(repro_counts_zips))) {
        for (i in seq_along(repro_counts_zips)) {
            temp_rds <- tempfile("repro")
            dir.create(temp_rds)
            unzipped = unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
            counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
            if (is.na(counts_file)) stop(paste("No *_counts.rds file found for index", i))
            counts_reprocessed <- as.data.frame(readRDS(counts_file))

            # ----- rdp16 -----
            if (!file.exists(file.path(local,"rdp16classified.csv.zip"))) {
            if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
                required_pkgs <- c("dada2", "Biostrings")
                missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
                if (length(missing_pkgs) > 0) {
                    stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                        ". Please install them before running this function.")
                }
                seqs <- Biostrings::DNAStringSet(colnames(counts_reprocessed))
                rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
                tax_reprocessed2 = make_taxa_label(rdpclassified) 
                tax_reprocessed2$Sequence <- sub("\\.\\.\\.[0-9]+$", "", rownames(tax_reprocessed2))
                rownames(tax_reprocessed2) <- tax_reprocessed2$Sequence
                tax_reprocessed2_list[[i]] <- tax_reprocessed2
                } else {
                stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
            }
            
            } else {
                tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
            }

            # Unzip and read taxonomy
            unzipped = unzip(repro_tax_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
            tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
            if (is.na(tax_file)) stop(paste("No *_taxa.rds file found for index", i))
            tax_reprocessed <- as.data.frame(readRDS(tax_file))
            tax_reprocessed <- make_taxa_label(tax_reprocessed)


            if (!raw) {
                aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = has_qPCR, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
                counts_reprocessed <- aligned$reprocessed
                counts_reprocessed2 = aligned$reprocessed
                matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
                matched_taxa2 <- tax_reprocessed2$Taxa[match(colnames(counts_reprocessed2), rownames(tax_reprocessed2))]
                colnames(counts_reprocessed) <- matched_taxa
                colnames(counts_reprocessed2) <- matched_taxa2
                counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
                counts_reprocessed2 <- collapse_duplicate_columns_exact(counts_reprocessed2)
                original_names <- colnames(counts_reprocessed)
                original_names2 <- colnames(counts_reprocessed2)
                counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
                counts_reprocessed2 <- as.data.frame(lapply(counts_reprocessed2, as.numeric), row.names = rownames(counts_reprocessed2), col.names = original_names2, check.names = FALSE)
                proportions_reprocessed2 <- sweep(counts_reprocessed2, 1, rowSums(counts_reprocessed2), '/')
            }
            # proportions
            proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

            # Label with study name based on zip filename prefix
            study_id <- sub("_.*$", "", basename(tools::file_path_sans_ext(repro_counts_zips[i])))
            tax_reprocessed$Study <- study_id
            tax_reprocessed2$Study <- study_id

            counts_reprocessed_list[[i]] <- counts_reprocessed
            proportions_reprocessed_list[[i]] <- proportions_reprocessed
            tax_reprocessed_list[[i]] <- tax_reprocessed

            counts_reprocessed2_list[[i]] <- counts_reprocessed2
            proportions_reprocessed2_list[[i]] <- proportions_reprocessed2

            cleanup_tempfiles(temp_rds)
        }
    

        # Combine all
        counts_reprocessed <- bind_rows(counts_reprocessed_list)
        proportions_reprocessed <- bind_rows(proportions_reprocessed_list)
        tax_reprocessed <- bind_rows(tax_reprocessed_list)

        if (!file.exists(file.path(local, "rdp16classified.csv.zip"))) {
            tax_reprocessed2 <- tax_reprocessed2_list %>%
                bind_rows() %>%
                distinct(Sequence, .keep_all = TRUE)
            proportions_reprocessed2 <- bind_rows(proportions_reprocessed2_list)
            counts_reprocessed2 <- bind_rows(counts_reprocessed2_list)
            write.csv(tax_reprocessed2, file = file.path(local, "rdp16classified.csv"), row.names = TRUE)
        } else {
            tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
            proportions_reprocessed2 <- bind_rows(proportions_reprocessed2_list) %>% as.data.frame()
            counts_reprocessed2 <- bind_rows(counts_reprocessed2_list) %>% as.data.frame()
        }
    }
    

    if (!raw) {
        counts_original = fill_na_zero_numeric(counts_original)
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_original = fill_na_zero_numeric(proportions_original)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
    }

    if (!raw) {
        metadata <- metadata %>%
            mutate(across(everything(), ~ifelse(. == "" | . == " ", NA, .))) %>%
            mutate(
                PatientID = factor(PatientID),
                Consistency = factor(replace_na(Consistency, "missing")),
                Pool = as.numeric(Pool),
                VanA = factor(replace_na(as.character(VanA), "missing"), levels = c("0", "1", "missing")),
                Accession = as.character(Accession),
                BioProject = as.character(BioProject),
                Timepoint = as.integer(Timepoint),
                DayRelativeToNearestHCT = as.integer(DayRelativeToNearestHCT)
            )
    }
    
    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original,
            reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
        ),
        tax = list(
            original = tax_original,
            reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
        ),
        scale = scale,
        metadata = metadata
    ))
}
