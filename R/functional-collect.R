#' Collect file-backed PICRUSt2 ASV contributions
#'
#' Reads the requested portion of a stratified PICRUSt2 table. Published tables
#' may be represented by one file or multiple shards; callers receive the same
#' data-frame schema in either case.
#'
#' @param x A PICRUSt2 branch returned inside a MUTT study result.
#' @param type Contribution type: `"ec"`, `"ko"`, or `"metacyc_abundance"`.
#' @param samples Optional exact sample identifiers.
#' @param functions Optional exact EC, KO, or MetaCyc identifiers.
#' @param asvs Optional original ASV identifiers from the source count table.
#' @param n_max Maximum number of matching rows returned.
#' @param collect_all Set `TRUE` to explicitly collect a complete sharded table.
#' @return A data frame containing PICRUSt2 contributions, original ASV
#' identifiers, and matched taxonomy.
#' @keywords internal
collect_functional <- function(
    x,
    type = c("ec", "ko", "metacyc_abundance"),
    samples = NULL,
    functions = NULL,
    asvs = NULL,
    n_max = Inf,
    collect_all = FALSE
) {
  read_picrust2_contributions(
    x,
    type = match.arg(type),
    samples = samples,
    functions = functions,
    asvs = asvs,
    n_max = n_max,
    collect_all = collect_all
  )
}

#' Convert file-backed PICRUSt2 contributions to a data frame
#'
#' @inheritParams collect_functional
#' @param row.names,optional Standard [as.data.frame()] arguments; unused.
#' @param ... Unused.
#' @export
as.data.frame.mutt_picrust_branch <- function(
    x,
    row.names = NULL,
    optional = FALSE,
    type = c("ec", "ko", "metacyc_abundance"),
    samples = NULL,
    functions = NULL,
    asvs = NULL,
    n_max = Inf,
    collect_all = FALSE,
    ...
) {
  collect_functional(
    x,
    type = match.arg(type),
    samples = samples,
    functions = functions,
    asvs = asvs,
    n_max = n_max,
    collect_all = collect_all
  )
}
