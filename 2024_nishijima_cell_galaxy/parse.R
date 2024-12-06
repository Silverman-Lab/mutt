require(stringr)

parse_2024_nishijima_cell_galaxy <- function() {
  local <- "2024_nishijima_cell_galaxy/"
  ## note read.table can handle compressed files nicely by default.
  out <- list()

  dat <- read.table(paste0(local, "GALAXY_load.tsv.gz"), header=TRUE)
  scale <- data.frame(count = dat$count)
  rownames(scale) <- dat$ID
  scale <- as(scale, "matrix")

  metadata <- data.frame(cohort=dat$cohort)
  rownames(metadata) <- dat$ID


  ## skipping first line because its terribly formated, will manually parse
  dat <- read.table(paste0(local, "GALAXY_mOTUs_v25.tsv.gz"),
                    header=FALSE, skip=1)
  proportions <- dat
  rownames(proportions) <- dat$V1
  proportions <- proportions[,-1]
  line1 <- readLines(paste0(local, "GALAXY_mOTUs_v25.tsv.gz"), n=1)
  taxnames <- strsplit(line1, "\t")[[1]]
  ## NOTE the taxa names in proportions are different than in the taxonomy table
  ## provided by authors. I am going with the taxonomy table provided by authors
  ## but this one looks like it could be useful as well.
  seqids <- str_extract(taxnames[-1], "(?<=mOTU_v25_)[:digit:]+")
  proportions <- as(t(proportions), "matrix")
  rownames(proportions) <- seqids
  ## TODO something is wrong with the second taxa in the table which is given a
  ## label of -1... I asked about it
  ## herehttps://github.com/grp-bork/microbial_load_predictor/issues/2

  ## TODO, still need to parse the tax table. But I ran out of time. 
  ## tax <- read.table(paste0(local, "motus2GTDB.txt.gz"))

  ## check sample names match
  all.equal(colnames(proportions), rownames(metadata))
  ## TODO more checking

  return(list(scale=scale, metadata=metadata, proportions=proportions))
 }
