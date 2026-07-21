.functional_publication_version <- 1L
.functional_default_max_asset_bytes <- 1900000000
.functional_default_shard_bytes <- 256 * 1024^2

.functional_validate_bytes <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1024) {
    stop("`", name, "` must be one finite number of at least 1024 bytes.", call. = FALSE)
  }
  as.numeric(value)
}

.functional_sha256 <- function(path) {
  unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

.functional_relative <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (!startsWith(path, prefix)) {
    stop("Published asset escaped its branch directory: ", path, call. = FALSE)
  }
  substring(path, nchar(prefix) + 1L)
}

.functional_asset_row <- function(path, root, role, source = "", rows = NA_real_) {
  data.frame(
    relative_path = .functional_relative(path, root),
    role = role,
    source = source,
    size_bytes = unname(file.info(path)$size),
    sha256 = .functional_sha256(path),
    rows = as.numeric(rows),
    stringsAsFactors = FALSE
  )
}

.functional_copy_asset <- function(source, destination, root, role, max_asset_bytes) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not publish functional asset: ", source, call. = FALSE)
  }
  size <- unname(file.info(destination)$size)
  if (!is.finite(size) || size > max_asset_bytes) {
    unlink(destination)
    stop(
      "Functional asset exceeds the per-file publication limit and is not a shardable table: ",
      source,
      call. = FALSE
    )
  }
  .functional_asset_row(destination, root, role, source)
}

.functional_open_text <- function(path, open = "rt") {
  if (grepl("[.]gz$", path, ignore.case = TRUE)) gzfile(path, open) else file(path, open)
}

.functional_line_values <- function(lines, column) {
  if (!length(lines) || is.na(column)) return(character())
  fields <- strsplit(lines, "\t", fixed = TRUE)
  values <- vapply(
    fields,
    function(x) if (length(x) >= column) x[[column]] else NA_character_,
    character(1)
  )
  sort(unique(values[!is.na(values) & nzchar(values)]))
}

.functional_write_tsv_piece <- function(lines, header, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(destination, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(header, lines), con, useBytes = TRUE)
  invisible(destination)
}

.functional_shard_tsv <- function(source, destination_dir, branch_root,
                                  max_asset_bytes, target_shard_bytes) {
  input <- .functional_open_text(source, "rt")
  on.exit(close(input), add = TRUE)
  header <- readLines(input, n = 1L, warn = FALSE)
  if (!length(header)) stop("Cannot publish an empty functional table: ", source, call. = FALSE)

  columns <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  sample_column <- match("sample", columns)
  function_column <- match("function", columns)
  taxon_column <- match("taxon", columns)
  if (is.na(sample_column)) sample_column <- NA_integer_
  if (is.na(function_column)) function_column <- NA_integer_
  if (is.na(taxon_column)) taxon_column <- NA_integer_

  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  shard_number <- 0L
  shard_rows <- list()
  asset_rows <- list()

  write_piece <- function(lines) {
    if (!length(lines)) return(invisible(NULL))
    shard_number <<- shard_number + 1L
    path <- file.path(destination_dir, sprintf("part-%05d.tsv.gz", shard_number))
    .functional_write_tsv_piece(lines, header, path)
    size <- unname(file.info(path)$size)
    if (size > max_asset_bytes) {
      unlink(path)
      shard_number <<- shard_number - 1L
      if (length(lines) == 1L) {
        stop("One functional table row exceeds the per-file publication limit.", call. = FALSE)
      }
      midpoint <- floor(length(lines) / 2L)
      write_piece(lines[seq_len(midpoint)])
      write_piece(lines[(midpoint + 1L):length(lines)])
      return(invisible(NULL))
    }

    samples <- .functional_line_values(lines, sample_column)
    functions <- .functional_line_values(lines, function_column)
    taxa <- .functional_line_values(lines, taxon_column)
    shard_rows[[length(shard_rows) + 1L]] <<- data.frame(
      relative_path = .functional_relative(path, branch_root),
      size_bytes = size,
      sha256 = .functional_sha256(path),
      rows = length(lines),
      stringsAsFactors = FALSE
    )
    shard_rows[[length(shard_rows)]]$sample_ids <- I(list(samples))
    shard_rows[[length(shard_rows)]]$function_ids <- I(list(functions))
    shard_rows[[length(shard_rows)]]$taxon_ids <- I(list(taxa))
    asset_rows[[length(asset_rows) + 1L]] <<- .functional_asset_row(
      path, branch_root, "stratified_shard", source, length(lines)
    )
    invisible(NULL)
  }

  target <- min(target_shard_bytes, max_asset_bytes * 0.75)
  buffer <- character()
  buffer_bytes <- 0
  repeat {
    lines <- readLines(input, n = 5000L, warn = FALSE)
    if (!length(lines)) break
    bytes <- nchar(lines, type = "bytes") + 1
    position <- 1L
    while (position <= length(lines)) {
      available <- target - buffer_bytes
      cumulative <- cumsum(bytes[position:length(lines)])
      take <- max(c(0L, which(cumulative <= available)))
      if (take == 0L && !length(buffer)) take <- 1L
      if (take > 0L) {
        selected <- position:(position + take - 1L)
        buffer <- c(buffer, lines[selected])
        buffer_bytes <- buffer_bytes + sum(bytes[selected])
        position <- position + take
      }
      if (buffer_bytes >= target || take == 0L) {
        write_piece(buffer)
        buffer <- character()
        buffer_bytes <- 0
      }
    }
  }
  write_piece(buffer)

  shards <- if (length(shard_rows)) do.call(rbind, shard_rows) else data.frame()
  assets <- if (length(asset_rows)) do.call(rbind, asset_rows) else data.frame()
  list(
    descriptor = list(
      storage = "sharded_tsv_gz",
      columns = columns,
      source_size_bytes = unname(file.info(source)$size),
      source_sha256 = .functional_sha256(source),
      rows = if (nrow(shards)) sum(shards$rows) else 0,
      shards = shards
    ),
    assets = assets
  )
}

.functional_publish_tsv <- function(source, destination_dir, branch_root,
                                    max_asset_bytes, target_shard_bytes,
                                    role = "functional_table") {
  if (!file.exists(source)) stop("Functional table is missing: ", source, call. = FALSE)
  source_size <- unname(file.info(source)$size)
  if (source_size <= max_asset_bytes) {
    destination <- file.path(destination_dir, basename(source))
    asset <- .functional_copy_asset(
      source, destination, branch_root, role, max_asset_bytes
    )
    con <- .functional_open_text(source, "rt")
    on.exit(close(con), add = TRUE)
    header <- readLines(con, n = 1L, warn = FALSE)
    columns <- if (length(header)) strsplit(header, "\t", fixed = TRUE)[[1L]] else character()
    return(list(
      descriptor = list(
        storage = "single_tsv_gz",
        relative_path = .functional_relative(destination, branch_root),
        columns = columns,
        size_bytes = source_size,
        sha256 = .functional_sha256(destination)
      ),
      assets = asset
    ))
  }
  .functional_shard_tsv(
    source, destination_dir, branch_root, max_asset_bytes, target_shard_bytes
  )
}

.functional_write_rds_component <- function(value, name, directory, branch_root,
                                            max_asset_bytes) {
  safe <- gsub("[^A-Za-z0-9._-]+", "_", name)
  path <- file.path(directory, paste0(safe, ".rds"))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path, version = 3)
  if (unname(file.info(path)$size) <= max_asset_bytes) {
    return(list(
      descriptor = list(storage = "single_rds", paths = .functional_relative(path, branch_root)),
      assets = .functional_asset_row(path, branch_root, "result_component", name)
    ))
  }
  unlink(path)
  if (!(is.matrix(value) || is.data.frame(value)) || nrow(value) < 2L) {
    stop("Result component cannot be safely divided below the publication limit: ", name, call. = FALSE)
  }

  pieces <- 2L
  repeat {
    groups <- split(seq_len(nrow(value)), cut(seq_len(nrow(value)), pieces, labels = FALSE))
    paths <- character(length(groups))
    sizes <- numeric(length(groups))
    for (index in seq_along(groups)) {
      paths[[index]] <- file.path(directory, sprintf("%s-part-%05d.rds", safe, index))
      saveRDS(value[groups[[index]], , drop = FALSE], paths[[index]], version = 3)
      sizes[[index]] <- unname(file.info(paths[[index]])$size)
    }
    if (all(sizes <= max_asset_bytes)) break
    unlink(paths)
    if (pieces >= nrow(value)) {
      stop("One row of result component exceeds the publication limit: ", name, call. = FALSE)
    }
    pieces <- min(nrow(value), pieces * 2L)
  }
  assets <- do.call(
    rbind,
    lapply(paths, .functional_asset_row, root = branch_root,
           role = "result_component_shard", source = name)
  )
  list(
    descriptor = list(
      storage = "row_sharded_rds",
      paths = vapply(paths, .functional_relative, character(1), root = branch_root),
      class = class(value),
      rows = nrow(value),
      columns = ncol(value)
    ),
    assets = assets
  )
}

.functional_write_result <- function(result, branch_root, max_asset_bytes) {
  path <- file.path(branch_root, "result.rds")
  saveRDS(result, path, version = 3)
  if (unname(file.info(path)$size) <= max_asset_bytes) {
    return(list(
      index = list(mode = "single", path = "result.rds"),
      assets = .functional_asset_row(path, branch_root, "result")
    ))
  }
  unlink(path)

  tabular <- names(result)[vapply(result, function(x) is.matrix(x) || is.data.frame(x), logical(1))]
  if (!length(tabular)) {
    stop("PICRUSt2 result.rds exceeds the publication limit and has no shardable top-level tables.", call. = FALSE)
  }
  skeleton <- result
  components <- list()
  assets <- list()
  for (name in tabular) {
    written <- .functional_write_rds_component(
      result[[name]], name, file.path(branch_root, "result_components"),
      branch_root, max_asset_bytes
    )
    components[[name]] <- written$descriptor
    assets[[length(assets) + 1L]] <- written$assets
    skeleton[[name]] <- NULL
  }
  skeleton_path <- file.path(branch_root, "result_skeleton.rds")
  saveRDS(skeleton, skeleton_path, version = 3)
  if (unname(file.info(skeleton_path)$size) > max_asset_bytes) {
    stop("Non-tabular PICRUSt2 result metadata exceeds the publication limit.", call. = FALSE)
  }
  assets[[length(assets) + 1L]] <- .functional_asset_row(
    skeleton_path, branch_root, "result_skeleton"
  )
  list(
    index = list(
      mode = "components",
      skeleton_path = "result_skeleton.rds",
      components = components,
      component_order = names(result)
    ),
    assets = do.call(rbind, assets)
  )
}

.functional_read_result <- function(branch_root, index) {
  if (identical(index$mode, "single")) {
    return(readRDS(file.path(branch_root, index$path)))
  }
  if (!identical(index$mode, "components")) {
    stop("Unknown functional result storage mode: ", index$mode, call. = FALSE)
  }
  result <- readRDS(file.path(branch_root, index$skeleton_path))
  for (name in names(index$components)) {
    component <- index$components[[name]]
    paths <- file.path(branch_root, unlist(component$paths, use.names = FALSE))
    values <- lapply(paths, readRDS)
    value <- if (identical(component$storage, "single_rds")) {
      values[[1L]]
    } else if (identical(component$storage, "row_sharded_rds")) {
      do.call(rbind, values)
    } else {
      stop("Unknown result component storage mode: ", component$storage, call. = FALSE)
    }
    result[[name]] <- value
  }
  result[index$component_order]
}

.functional_is_shardable_table <- function(path) {
  grepl("[.](tsv|txt)([.]gz)?$", path, ignore.case = TRUE)
}

.functional_publish_branch <- function(source_root, branch_source, branch_target,
                                       max_asset_bytes, target_shard_bytes) {
  relative_branch <- .functional_relative(branch_source, source_root)
  parts <- strsplit(relative_branch, "/", fixed = TRUE)[[1L]]
  method <- parts[[1L]]
  branch <- paste(parts[-1L], collapse = "/")
  result <- readRDS(file.path(branch_source, "result.rds"))
  dir.create(branch_target, recursive = TRUE, showWarnings = FALSE)
  assets <- list()
  contribution_sources <- character()

  if (identical(method, "picrust2")) {
    if (!is.list(result$stratified) || !length(result$stratified)) {
      stop("PICRUSt2 result has no stratified descriptors: ", relative_branch, call. = FALSE)
    }
    for (type in names(result$stratified)) {
      descriptor <- result$stratified[[type]]
      source <- file.path(branch_source, descriptor$relative_path)
      contribution_sources <- c(contribution_sources, normalizePath(source, winslash = "/", mustWork = TRUE))
      published <- .functional_publish_tsv(
        source,
        file.path(branch_target, "stratified", .safe_branch_name(type)),
        branch_target,
        max_asset_bytes,
        target_shard_bytes,
        role = paste0("stratified_", type)
      )
      if (!is.null(descriptor$columns)) {
        published$descriptor$columns <- descriptor$columns
      }
      result$stratified[[type]] <- published$descriptor
      assets[[length(assets) + 1L]] <- published$assets
    }
  }

  source_files <- list.files(
    branch_source, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE
  )
  source_files <- source_files[file.info(source_files)$isdir %in% FALSE]
  normalized <- normalizePath(source_files, winslash = "/", mustWork = TRUE)
  relative_files <- substring(
    normalized,
    nchar(normalizePath(branch_source, winslash = "/", mustWork = TRUE)) + 2L
  )
  keep <- basename(relative_files) != "result.rds" &
    !grepl("(^|/)raw/intermediate(/|$)", relative_files) &
    !grepl("(^|/)input(/|$)", relative_files) &
    !(normalized %in% contribution_sources)

  for (index in which(keep)) {
    source <- source_files[[index]]
    relative <- relative_files[[index]]
    destination <- file.path(branch_target, relative)
    if (unname(file.info(source)$size) <= max_asset_bytes) {
      assets[[length(assets) + 1L]] <- .functional_copy_asset(
        source, destination, branch_target, "retained_output", max_asset_bytes
      )
    } else if (.functional_is_shardable_table(source)) {
      published <- .functional_publish_tsv(
        source,
        file.path(branch_target, "oversized", .safe_branch_name(relative)),
        branch_target,
        max_asset_bytes,
        target_shard_bytes,
        role = "retained_output_shard"
      )
      assets[[length(assets) + 1L]] <- published$assets
    } else {
      stop("Oversized non-tabular functional output cannot be published safely: ", source, call. = FALSE)
    }
  }

  result$provenance$output_directory <- branch_target
  written_result <- .functional_write_result(result, branch_target, max_asset_bytes)
  assets[[length(assets) + 1L]] <- written_result$assets
  list(
    branch = list(
      method = method,
      branch = branch,
      relative_directory = relative_branch,
      result = written_result$index
    ),
    assets = do.call(rbind, assets)
  )
}

#' Build a size-safe functional-data publication
#'
#' Constructs the API-facing `functional_data/` directory from a validated
#' execution cache without modifying the cache. Oversized tabular outputs are
#' deterministically divided into gzip-compressed shards.
#'
#' @param study_dir One MUTT study directory containing `functional/`.
#' @param max_asset_bytes Maximum allowed size of one published file.
#' @param target_shard_bytes Target uncompressed bytes per tabular shard.
#' @return Path to the published `functional_data/` directory, invisibly.
#' @keywords internal
publish_functional_data <- function(
    study_dir,
    max_asset_bytes = .functional_default_max_asset_bytes,
    target_shard_bytes = min(.functional_default_shard_bytes, max_asset_bytes * 0.75)
) {
  max_asset_bytes <- .functional_validate_bytes(max_asset_bytes, "max_asset_bytes")
  target_shard_bytes <- .functional_validate_bytes(target_shard_bytes, "target_shard_bytes")
  source_root <- file.path(study_dir, "functional")
  if (!dir.exists(source_root)) {
    stop("Validated functional execution cache not found: ", source_root, call. = FALSE)
  }

  target <- file.path(study_dir, "functional_data")
  staging <- tempfile("functional-data-staging-", tmpdir = study_dir)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)

  result_paths <- list.files(
    source_root, pattern = "^result[.]rds$", recursive = TRUE, full.names = TRUE
  )
  result_paths <- result_paths[
    grepl("(^|/)(picrust2|faprotax)/", gsub("\\\\", "/", result_paths))
  ]
  if (!length(result_paths)) {
    stop("No validated PICRUSt2 or FAPROTAX result.rds files were found.", call. = FALSE)
  }

  branches <- list()
  assets <- list()
  for (path in sort(result_paths)) {
    branch_source <- dirname(path)
    relative <- .functional_relative(branch_source, source_root)
    branch_target <- file.path(staging, relative)
    published <- .functional_publish_branch(
      source_root, branch_source, branch_target, max_asset_bytes, target_shard_bytes
    )
    branches[[length(branches) + 1L]] <- published$branch
    branch_assets <- published$assets
    branch_assets$branch <- relative
    assets[[length(assets) + 1L]] <- branch_assets
  }
  assets <- do.call(rbind, assets)

  functional_manifest <- file.path(source_root, "manifest.json")
  if (file.exists(functional_manifest)) {
    file.copy(functional_manifest, file.path(staging, "functional_manifest.json"), overwrite = TRUE)
  }
  index <- list(
    schema_version = .functional_publication_version,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    branches = branches,
    assets = assets
  )
  saveRDS(index, file.path(staging, "publication_index.rds"), version = 3)
  if (unname(file.info(file.path(staging, "publication_index.rds"))$size) > max_asset_bytes) {
    stop("Functional publication index exceeds the per-file limit.", call. = FALSE)
  }
  jsonlite::write_json(
    list(
      schema_version = index$schema_version,
      created_utc = index$created_utc,
      branches = lapply(branches, function(x) x[c("method", "branch", "relative_directory")]),
      assets = assets
    ),
    file.path(staging, "publication_manifest.json"),
    dataframe = "rows", auto_unbox = TRUE, pretty = TRUE, na = "null"
  )

  published_files <- list.files(staging, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  published_files <- published_files[file.info(published_files)$isdir %in% FALSE]
  oversized <- published_files[unname(file.info(published_files)$size) > max_asset_bytes]
  if (length(oversized)) {
    stop("Published functional assets exceed the per-file limit: ", paste(oversized, collapse = ", "), call. = FALSE)
  }

  backup <- paste0(target, ".previous")
  if (dir.exists(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  if (dir.exists(target) && !file.rename(target, backup)) {
    stop("Could not protect the previous functional_data publication.", call. = FALSE)
  }
  if (!file.rename(staging, target)) {
    if (dir.exists(backup)) file.rename(backup, target)
    stop("Could not publish functional_data atomically.", call. = FALSE)
  }
  if (dir.exists(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  invisible(target)
}

.functional_validate_publication_assets <- function(root, index) {
  assets <- index$assets
  required <- c("branch", "relative_path", "size_bytes")
  if (!is.data.frame(assets) || !all(required %in% names(assets))) {
    stop("Functional publication index has an invalid asset table.", call. = FALSE)
  }
  paths <- file.path(
    root,
    as.character(assets$branch),
    as.character(assets$relative_path)
  )
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop(
      "Functional publication is incomplete; missing asset(s): ",
      paste(paths[missing], collapse = ", "),
      call. = FALSE
    )
  }
  observed <- unname(file.info(paths)$size)
  expected <- as.numeric(assets$size_bytes)
  changed <- is.na(observed) | is.na(expected) | observed != expected
  if (any(changed)) {
    stop(
      "Functional publication asset size does not match its manifest: ",
      paste(paths[changed], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.load_functional_publication <- function(study_dir) {
  root <- file.path(study_dir, "functional_data")
  index_path <- file.path(root, "publication_index.rds")
  if (!file.exists(index_path)) return(NULL)
  index <- readRDS(index_path)
  if (!identical(as.integer(index$schema_version), .functional_publication_version)) {
    stop("Unsupported functional publication schema version.", call. = FALSE)
  }
  .functional_validate_publication_assets(root, index)
  out <- .empty_functional_result(TRUE)
  for (entry in index$branches) {
    branch_root <- file.path(root, entry$relative_directory)
    result <- .functional_read_result(branch_root, entry$result)
    if (is.list(result$provenance)) result$provenance$output_directory <- branch_root
    path <- strsplit(entry$branch, "/", fixed = TRUE)[[1L]]
    if (identical(entry$method, "picrust2")) {
      class(result) <- unique(c("mutt_picrust_branch", class(result)))
      out$picrust2 <- .set_nested(out$picrust2, path, result)
    } else if (identical(entry$method, "faprotax")) {
      out$faprotax <- .set_nested(out$faprotax, path, result)
    }
  }
  manifest_path <- file.path(root, "functional_manifest.json")
  if (file.exists(manifest_path)) {
    manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
    out$manifest <- as.data.frame(manifest, stringsAsFactors = FALSE)
  }
  attr(out, "publication") <- index[c("schema_version", "created_utc")]
  out
}

.functional_descriptor_shards <- function(descriptor) {
  storage <- if (is.null(descriptor$storage)) "legacy_single" else descriptor$storage
  if (storage %in% c("single_tsv_gz", "legacy_single")) {
    return(data.frame(
      relative_path = descriptor$relative_path,
      stringsAsFactors = FALSE
    ))
  }
  if (!identical(storage, "sharded_tsv_gz") || !is.data.frame(descriptor$shards)) {
    stop("Unknown PICRUSt2 contribution storage descriptor.", call. = FALSE)
  }
  descriptor$shards
}

.functional_filter_shards <- function(shards, samples = NULL,
                                      functions = NULL, taxa = NULL) {
  if (!nrow(shards)) return(shards)
  keep <- rep(TRUE, nrow(shards))
  selectors <- list(
    sample_ids = samples,
    function_ids = functions,
    taxon_ids = taxa
  )
  for (column in names(selectors)) {
    selected <- selectors[[column]]
    if (is.null(selected) || !column %in% names(shards)) next
    keep <- keep & vapply(
      shards[[column]],
      function(ids) any(selected %in% ids),
      logical(1)
    )
  }
  shards[keep, , drop = FALSE]
}

.read_functional_descriptor <- function(branch_root, descriptor, samples = NULL,
                                        functions = NULL, taxa = NULL,
                                        n_max = Inf, collect_all = FALSE) {
  storage <- if (is.null(descriptor$storage)) "legacy_single" else descriptor$storage
  if (identical(storage, "sharded_tsv_gz") && is.infinite(n_max) &&
      is.null(samples) && is.null(functions) && is.null(taxa) && !isTRUE(collect_all)) {
    stop(
      "This contribution table is sharded. Supply `samples`, `functions`, `asvs`, ",
      "or set `collect_all = TRUE` to read the complete table.",
      call. = FALSE
    )
  }
  shards <- .functional_filter_shards(
    .functional_descriptor_shards(descriptor),
    samples = samples,
    functions = functions,
    taxa = taxa
  )
  if (!nrow(shards)) {
    return(as.data.frame(matrix(nrow = 0L, ncol = length(descriptor$columns),
                                dimnames = list(NULL, descriptor$columns))))
  }
  remaining <- n_max
  output <- list()
  for (index in seq_len(nrow(shards))) {
    path <- file.path(branch_root, shards$relative_path[[index]])
    if (!file.exists(path)) stop("Published PICRUSt2 shard is missing: ", path, call. = FALSE)
    filtered <- !is.null(samples) || !is.null(functions) || !is.null(taxa)
    read_limit <- if (filtered || is.infinite(remaining)) Inf else as.integer(remaining)
    table <- readr::read_tsv(
      path, n_max = read_limit, show_col_types = FALSE,
      progress = FALSE, name_repair = "minimal"
    )
    if (!is.null(samples)) table <- table[table$sample %in% samples, , drop = FALSE]
    if (!is.null(functions)) table <- table[table[["function"]] %in% functions, , drop = FALSE]
    if (!is.null(taxa)) table <- table[table$taxon %in% taxa, , drop = FALSE]
    if (!is.infinite(remaining) && nrow(table) > remaining) {
      table <- table[seq_len(remaining), , drop = FALSE]
    }
    output[[length(output) + 1L]] <- table
    if (!is.infinite(remaining)) {
      remaining <- remaining - nrow(table)
      if (remaining <= 0) break
    }
  }
  if (!length(output)) {
    return(as.data.frame(matrix(nrow = 0L, ncol = length(descriptor$columns),
                                dimnames = list(NULL, descriptor$columns))))
  }
  as.data.frame(do.call(rbind, output), check.names = FALSE, stringsAsFactors = FALSE)
}

read_picrust2_contributions <- function(
    x,
    type = c("ec", "ko", "metacyc_abundance"),
    samples = NULL,
    functions = NULL,
    asvs = NULL,
    n_max = Inf,
    collect_all = FALSE
) {
  type <- match.arg(type)
  if (!is.list(x) || is.null(x$stratified) || is.null(x$asv_mapping) ||
      is.null(x$provenance$output_directory)) {
    stop("`x` must be a PICRUSt2 branch result returned by MUTT.", call. = FALSE)
  }
  if (!is.numeric(n_max) || length(n_max) != 1L || is.na(n_max) || n_max <= 0) {
    stop("`n_max` must be a positive number or Inf.", call. = FALSE)
  }
  for (argument in c("samples", "functions", "asvs")) {
    value <- get(argument)
    if (!is.null(value) && (!is.character(value) || anyNA(value))) {
      stop("`", argument, "` must be NULL or a character vector.", call. = FALSE)
    }
  }
  if (!is.logical(collect_all) || length(collect_all) != 1L || is.na(collect_all)) {
    stop("`collect_all` must be TRUE or FALSE.", call. = FALSE)
  }

  picrust_taxa <- NULL
  if (!is.null(asvs)) {
    picrust_taxa <- x$asv_mapping$picrust_id[
      match(asvs, x$asv_mapping$original_feature_id)
    ]
    if (anyNA(picrust_taxa)) {
      missing <- asvs[is.na(picrust_taxa)]
      stop("Unknown original ASV identifier(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }
  descriptor <- x$stratified[[type]]
  out <- .read_functional_descriptor(
    x$provenance$output_directory,
    descriptor,
    samples = samples,
    functions = functions,
    taxa = picrust_taxa,
    n_max = n_max,
    collect_all = collect_all
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
  cbind(out, taxonomy)
}
