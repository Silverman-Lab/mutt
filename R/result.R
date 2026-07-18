.new_mutt_result <- function(data, audit, provenance, validation) {
  class(data) <- c("mutt_result", "list")
  attr(data, "audit") <- audit
  attr(data, "provenance") <- provenance
  attr(data, "validation") <- validation
  data
}

#' Print a MUTT result summary
#'
#' @param x A result returned by [mutt()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mutt_result <- function(x, ...) {
  audit <- attr(x, "audit")
  cat("<mutt_result>\n")
  cat("  requested:", if (is.null(audit)) length(x) else nrow(audit), "\n")
  cat("  successful:", sum(!vapply(x, is.null, logical(1))), "\n")
  if (!is.null(audit) && any(audit$status != "success")) {
    cat("  failed:", sum(audit$status != "success"), "\n")
  }
  cat("  studies:", paste(names(x)[!vapply(x, is.null, logical(1))], collapse = ", "), "\n")
  invisible(x)
}

#' Read a file-backed PICRUSt2 contribution table
#'
#' @param x A PICRUSt2 branch returned inside a MUTT study result.
#' @param row.names,optional Standard [as.data.frame()] arguments; unused.
#' @param type Contribution type: `"ec"`, `"ko"`, or `"metacyc_abundance"`.
#' @param n_max Maximum rows to read.
#' @param ... Unused.
#' @return A data frame joined to original ASV identifiers and taxonomy.
#' @export
as.data.frame.mutt_picrust_branch <- function(x, row.names = NULL, optional = FALSE,
                                               type = c("ec", "ko", "metacyc_abundance"),
                                               n_max = Inf, ...) {
  read_picrust2_contributions(x, type = match.arg(type), n_max = n_max)
}
