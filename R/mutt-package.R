#' MUTT package internals
#'
#' @keywords internal
#' @import dplyr
#' @import optparse
#' @import stringr
#' @import tidyr
#' @import tidyverse
#' @importFrom stats setNames
#' @importFrom utils read.csv read.table unzip
"_PACKAGE"

utils::globalVariables(c(
  "Study Identifier", "OTU_ID", "qiime_sklearn", "taxonomy",
  "D0", "D1", "D2", "D3", "D4", "D5", "D6", "Kingdom", "Species",
  "usearch_sintax_80%", "usearch_sintax", "usearch_utax",
  "vsearch_usearchglobal"
))
