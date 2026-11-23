add_sequence_column <-
function (x) 
{
    if (is.data.frame(x)) {
        seq_pattern <- "^[ACGTRYSWKMBDHVN]+$"
        seq_rows <- grepl(seq_pattern, rownames(x), ignore.case = TRUE)
        if ("sequence" %in% names(x)) {
            x$Sequence <- as.character(x$sequence)
            x$sequence <- NULL
        }
        if (!"Sequence" %in% names(x)) {
            x$Sequence <- rep(NA_character_, nrow(x))
        }
        else {
            x$Sequence <- as.character(x$Sequence)
        }
        valid_sequence <- !is.na(x$Sequence) & grepl(seq_pattern, 
            x$Sequence, ignore.case = TRUE)
        x$Sequence[!valid_sequence] <- NA_character_
        to_fill <- seq_rows & is.na(x$Sequence)
        x$Sequence[to_fill] <- rownames(x)[to_fill]
        return(x)
    }
    if (is.list(x)) {
        return(lapply(x, add_sequence_column))
    }
    x
}
build_taxonomy_table <-
function (df, method = c("qiime", "sintax", "sintax_full", "utax", 
    "vsearch")) 
{
    method <- match.arg(method)
    if (method == "qiime") {
        tax_source <- df %>% select(OTU_ID, taxonomy = qiime_sklearn)
        if (any(grepl("^D_\\d+__", tax_source$taxonomy, perl = TRUE))) {
            tax_df <- tax_source %>% separate(taxonomy, into = c("D0", 
                "D1", "D2", "D3", "D4", "D5", "D6"), sep = ";", 
                fill = "right") %>% transmute(OTU_ID, Kingdom = gsub("^D_0__", 
                "", D0), Phylum = gsub("^D_1__", "", D1), Class = gsub("^D_2__", 
                "", D2), Order = gsub("^D_3__", "", D3), Family = gsub("^D_4__", 
                "", D4), Genus = gsub("^D_5__", "", D5), Species = gsub("^D_6__", 
                "", D6))
        }
        else {
            tax_df <- tax_source %>% separate(taxonomy, into = c("Kingdom", 
                "Phylum", "Class", "Order", "Family", "Genus", 
                "Species"), sep = ";", fill = "right") %>% mutate(across(Kingdom:Species, 
                ~gsub("^[a-z]__", "", .)))
        }
        tax_df <- tax_df %>% mutate(across(Kingdom:Species, ~ifelse(. == 
            "" | is.na(.), NA, .)))
    }
    else if (method == "sintax") {
        tax_df <- df %>% select(OTU_ID, taxonomy = `usearch_sintax_80%`) %>% 
            separate(taxonomy, into = c("Kingdom", "Phylum", 
                "Class", "Order", "Family", "Genus", "Species"), 
                sep = ",", fill = "right") %>% mutate(across(Kingdom:Species, 
            ~gsub("^[a-z]:", "", .)))
    }
    else if (method == "sintax_full") {
        tax_df <- df %>% select(OTU_ID, taxonomy = usearch_sintax) %>% 
            mutate(taxonomy = gsub("\\([^)]+\\)", "", taxonomy)) %>% 
            separate(taxonomy, into = c("Kingdom", "Phylum", 
                "Class", "Order", "Family", "Genus", "Species"), 
                sep = ",", fill = "right") %>% mutate(across(Kingdom:Species, 
            ~gsub("^[a-z]:", "", .)))
    }
    else if (method == "utax") {
        tax_df <- df %>% select(OTU_ID, taxonomy = usearch_utax) %>% 
            mutate(taxonomy = gsub(".*\\|refs\\|", "", taxonomy)) %>% 
            separate(taxonomy, into = c("Kingdom", "Phylum", 
                "Class", "Order", "Family", "Genus", "Species"), 
                sep = ";", fill = "right") %>% mutate(across(Kingdom:Species, 
            ~gsub("^[a-z]__", "", .)))
    }
    else if (method == "vsearch") {
        tax_df <- df %>% select(OTU_ID, taxonomy = vsearch_usearchglobal) %>% 
            mutate(taxonomy = gsub(".*\\|refs\\|", "", taxonomy)) %>% 
            separate(taxonomy, into = c("Kingdom", "Phylum", 
                "Class", "Order", "Family", "Genus", "Species"), 
                sep = ";", fill = "right") %>% mutate(across(Kingdom:Species, 
            ~gsub("^[a-z]__", "", .)))
    }
    tax_df <- tax_df %>% mutate(across(Kingdom:Species, ~ifelse(. == 
        "" | is.na(.), NA, .)))
    return(tax_df)
}
cleanup_tempfiles <-
function (temp_paths, safety_check = TRUE) 
{
    for (p in temp_paths) {
        if (!file.exists(p)) 
            next
        if (normalizePath(p) == normalizePath(tempdir())) {
            warning("Not deleting session-wide tempdir(): ", 
                p)
            next
        }
        if (dir.exists(p) && safety_check) {
            if (grepl("^(/tmp|/var/tmp|/private/tmp)", normalizePath(p))) {
                unlink(p, recursive = TRUE, force = TRUE)
            }
            else {
                warning("Skipping deletion of non-tmp directory: ", 
                  p)
            }
        }
        else {
            unlink(p, recursive = TRUE, force = TRUE)
        }
    }
}
collapse_duplicate_columns_exact <-
function (df) 
{
    orig_rn <- rownames(df)
    orig_cn <- colnames(df)
    orig_cn[is.na(orig_cn)] <- "Unclassified"
    df_num <- data.frame(lapply(df, function(col) {
        if (is.numeric(col)) {
            col
        }
        else {
            as.numeric(as.character(col))
        }
    }), row.names = orig_rn, check.names = FALSE)
    summed <- rowsum(t(as.matrix(df_num)), group = orig_cn, reorder = FALSE)
    result <- as.data.frame(t(summed), row.names = orig_rn, check.names = FALSE)
    return(result)
}
fill_na_zero_numeric <-
function (x) 
{
    if (missing(x)) 
        return(NULL)
    if (is.data.frame(x)) {
        x[] <- lapply(x, function(y) if (is.numeric(y)) 
            replace(y, is.na(y), 0)
        else y)
    }
    else if (is.matrix(x) && is.numeric(x)) {
        x[is.na(x)] <- 0
    }
    else if (is.list(x)) {
        x <- lapply(x, fill_na_zero_numeric)
    }
    x
}
first_row_to_colnames <-
function (df) 
{
    colnames(df) <- as.character(df[1, ])
    df <- df[-1, ]
    return(df)
}
make_taxa_label <-
function (df) 
{
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
        "Genus")
    prefixes <- c("k", "p", "c", "o", "f", "g")
    has_sp <- "Species" %in% colnames(df)
    if (has_sp) {
        tax_ranks <- c(tax_ranks, "Species")
        prefixes <- c(prefixes, "s")
    }
    if (!all(tax_ranks %in% colnames(df))) {
        stop("Data frame must contain columns: ", paste(tax_ranks, 
            collapse = ", "))
    }
    romans <- c("i", "ii", "iii", "iv", "v", "vi", "vii", "viii", 
        "ix", "x", "xi", "xii", "xiii", "xiv", "xv", "xvi", "xvii", 
        "xviii", "xix", "xx")
    df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
        x <- as.character(x)
        x <- gsub("^[a-z]_{1,2}", "", x)
        x[is.na(x) | trimws(x) == ""] <- "unclassified"
        for (rn in romans) {
            x <- gsub(paste0("(?i)\\b", rn, "\\b"), toupper(rn), 
                x, perl = TRUE)
        }
        x <- gsub("[[:space:][:punct:]]+", "_", x)
        x <- gsub("_+", "_", x)
        x <- gsub("^_|_$", "", x)
        x
    })
    df$Taxa <- apply(df[, tax_ranks], 1, function(tr) {
        if (has_sp && tr["Species"] != "unclassified") {
            return(paste0("s_", tr["Species"]))
        }
        if (tr["Genus"] != "unclassified") {
            return(paste0("g_", tr["Genus"]))
        }
        fallback_start <- length(prefixes) - ifelse(has_sp, 2, 
            1)
        for (i in seq(fallback_start, 1)) {
            val <- tr[tax_ranks[i]]
            if (val != "unclassified") {
                return(paste0("uc_", prefixes[i], "_", val))
            }
        }
        "unclassified"
    })
    df
}
read_xlsx_zip <-
function (zipfile, sheet = NULL, skip = 0, tmp = tempdir()) 
{
    if (!file.exists(zipfile)) {
        stop("Zip file not found: ", zipfile)
    }
    zinfo <- utils::unzip(zipfile, list = TRUE)
    xlsx_name <- zinfo$Name[grepl("\\.xlsx$", zinfo$Name)][1]
    if (is.na(xlsx_name)) {
        stop("No .xlsx file found inside ", zipfile)
    }
    extracted <- utils::unzip(zipfile, files = xlsx_name, exdir = tmp, 
        overwrite = TRUE, junkpaths = TRUE)
    xlsx_path <- file.path(tmp, basename(xlsx_name))
    if (!file.exists(xlsx_path)) {
        stop("Extraction failed--'", xlsx_path, "' does not exist")
    }
    dat <- readxl::read_xlsx(xlsx_path, sheet = sheet, skip = skip)
    return(dat)
}
read_zipped_csv <-
function (zip_path) 
{
    if (file.exists(zip_path)) {
        csv_file <- unzip(zip_path, list = TRUE)$Name[1]
        read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
    }
    else {
        warning(paste("File not found:", zip_path))
        return(NA)
    }
}
read_zipped_table <-
function (zip_path, sep = ",", header = TRUE, row.names = 1, 
    check.names = FALSE, skip = 0) 
{
    if (file.exists(zip_path)) {
        inner_file <- unzip(zip_path, list = TRUE)$Name[1]
        con <- unz(zip_path, inner_file)
        if (is.null(row.names)) {
            read.table(con, sep = sep, header = header, row.names = NULL, 
                check.names = check.names, stringsAsFactors = FALSE, 
                skip = skip)
        }
        else {
            read.table(con, sep = sep, header = header, row.names = row.names, 
                check.names = check.names, stringsAsFactors = FALSE, 
                skip = skip)
        }
    }
    else {
        warning(paste("File not found:", zip_path))
        return(NA)
    }
}
remove_empty_columns <-
function (df) 
{
    if (!is.data.frame(df)) {
        stop("Input must be a data frame")
    }
    is_empty_column <- function(x) {
        if (is.numeric(x)) {
            all(is.na(x))
        }
        else {
            all(is.na(x) | trimws(x) == "")
        }
    }
    df[, !sapply(df, is_empty_column), drop = FALSE]
}
rename_and_align <-
function (counts_reprocessed = NULL, counts_original = NULL, 
    proportions_original = NULL, metadata, scale, by_col, align = TRUE, 
    study_name = NULL) 
{
    stopifnot(is.character(by_col))
    if (!all(by_col %in% colnames(scale))) {
        message("\n", study_name, " | by_col '", by_col, "' not found in scale columns. Returning original dataframes.")
        return(list(reprocessed = counts_reprocessed, counts_original = counts_original, 
            proportions_original = proportions_original))
    }
    scale[by_col] <- lapply(scale[by_col], as.character)
    target_samples <- if (align) {
        unique(do.call(paste, c(scale[by_col], sep = "_")))
    }
    else NULL
    if (!is.null(counts_reprocessed) && (is.matrix(counts_reprocessed) || 
        is.data.frame(counts_reprocessed))) {
        stopifnot("Accession" %in% colnames(metadata), all(by_col %in% 
            colnames(metadata)))
        metadata[by_col] <- lapply(metadata[by_col], as.character)
        acc_to_name <- setNames(do.call(paste, c(metadata[by_col], 
            sep = "_")), metadata$Accession)
        if (length(intersect(rownames(counts_reprocessed), metadata$Accession)) == 
            0) {
            warning("\n", study_name, " | No matching accessions found between counts_reprocessed and metadata")
            reprocessed = counts_reprocessed
        }
        else {
            mapped_names <- acc_to_name[rownames(counts_reprocessed)]
            valid <- !is.na(mapped_names)
            dropped <- sum(!valid)
            if (dropped > 0) 
                message("\n", study_name, " | counts_reprocessed: dropped ", 
                  dropped, " unmatched accessions")
            counts_reprocessed <- counts_reprocessed[valid, , 
                drop = FALSE]
            mapped_names <- make.unique(mapped_names[valid])
            rownames(counts_reprocessed) <- mapped_names
        }
    }
    if (align) {
        if (!is.null(counts_reprocessed)) {
            before_n <- nrow(counts_reprocessed)
            counts_reprocessed <- counts_reprocessed[rownames(counts_reprocessed) %in% 
                target_samples, , drop = FALSE]
            after_n <- nrow(counts_reprocessed)
            if (before_n != after_n) 
                message("\n", study_name, " | counts_reprocessed: dropped ", 
                  before_n - after_n, " unaligned samples")
        }
        if (!is.null(counts_original) && (is.matrix(counts_original) || 
            is.data.frame(counts_original))) {
            before_n <- nrow(counts_original)
            counts_original <- counts_original[rownames(counts_original) %in% 
                target_samples, , drop = FALSE]
            after_n <- nrow(counts_original)
            if (before_n != after_n) 
                message("\n", study_name, " | counts_original: dropped ", 
                  before_n - after_n, " unaligned samples")
        }
        if (!is.null(proportions_original) && (is.matrix(proportions_original) || 
            is.data.frame(proportions_original))) {
            before_n <- nrow(proportions_original)
            proportions_original <- proportions_original[rownames(proportions_original) %in% 
                target_samples, , drop = FALSE]
            after_n <- nrow(proportions_original)
            if (before_n != after_n) 
                message("\n", study_name, " | proportions_original: dropped ", 
                  before_n - after_n, " unaligned samples")
        }
    }
    list(reprocessed = counts_reprocessed, counts_original = counts_original, 
        proportions_original = proportions_original)
}
standardize_output_order <-
function (output) 
{
    if (!is.list(output)) {
        stop("Input must be a list")
    }
    standard_order <- c("counts", "proportions", "tax", "function", 
        "scale", "metadata")
    current_elements <- names(output)
    other_elements <- setdiff(current_elements, standard_order)
    result <- list()
    for (elem in standard_order) {
        if (elem %in% current_elements) {
            result[[elem]] <- output[[elem]]
        }
    }
    for (elem in other_elements) {
        result[[elem]] <- output[[elem]]
    }
    result
}
# Final harmonization previously performed by mutt(), retained with the recipes.
remove_all_na_fields <- function(x) {
  if(is.data.frame(x)) return(x[,vapply(x,function(column)!all(is.na(column)),logical(1)),drop=FALSE])
  if(is.list(x)) return(lapply(x,remove_all_na_fields))
  x
}
