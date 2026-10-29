# mutt

`mutt` is a curated library of microbiome studies pairing taxonomic profiles with
microbial-load measurements. Install a small R package, list studies, and request
only the processed datasets you need.

```r
install.packages("remotes")
remotes::install_github("Silverman-Lab/mutt")
library(mutt)
mutt_studies()
x <- mutt("vandeputte2017")
y <- mutt(c("vandeputte2017", "contijoch2019", "morton2019"))
```

**Release preparation status:** this branch contains the new distribution API.
The bundled manifest explicitly distinguishes `staged` assets from public
`available` assets and studies with no usable current output (`unavailable`).
Staged data require the maintainer upload/publication step before normal downloads
work. `mutt()` explains this instead of following a fabricated or mutable URL.

Only requested studies are downloaded. You do not need a repository clone, Git
LFS, Docker, a bioinformatics environment or any preprocessing tools. Each study
is pinned to a versioned GitHub Release asset, byte size and SHA-256 checksum.
Verified data are cached under `tools::R_user_dir("mutt", "cache")`.

```r
x <- mutt("vandeputte2017")                         # reuses verified cache
x <- mutt("vandeputte2017", refresh = TRUE)         # re-downloads and verifies
mutt_cache_remove("vandeputte2017")                # explicit removal
```

`mutt()` always returns a named list, including for one study. A study preserves
its existing `counts`, `proportions`, `tax`, `scale` and `metadata` fields, other
scientific fields, and all named nested assay/classifier/cohort branches. Load
measurements retain their original column names and units. Some fields contain
NA or heterogeneous representations because the contributing studies differ.
Choose a representation explicitly for your analysis; no scientific branch is
silently selected or flattened.

The v1 baseline is the output of the previous package at commit `01c648b`, reading
its stored processed tables, with `rawdata=FALSE`, `align_samples=FALSE` and
`functional=FALSE`. The empty inferred-function field is intentionally removed.
Counts, taxon names, sample identifiers, scale, taxonomy and metadata remain
identical. The old sample-aligned output is stored separately:

```r
aligned <- mutt("vandeputte2017", representation = "aligned")
```

Previous parsers that already failed are listed with their failure reason in the
migration audit. An older saved object is never silently substituted for current
data. No DADA2 denoising, taxonomy classification or MetaPhlAn profiling was run
as part of migration.

## Sources and preprocessing

```r
mutt_source("vandeputte2017")                      # sources, notes, accessions
mutt_recipe("vandeputte2017")                      # pinned recipe URLs
mutt_download_raw("vandeputte2017", dest = "inputs") # stored source inputs
# Explicitly request public FASTQs too; this can be a large download:
# mutt_download_raw("vandeputte2017", dest = "inputs", reads = TRUE)
```

Source inputs include published tables, supplementary spreadsheets and existing
processed taxonomic inputs. Exact original inputs are available through
checksum-pinned historical archive URLs; authoritative references/accessions are
recorded beside them. Raw sequencing reads are retrieved from ENA when explicitly
requested, never mirrored in the data release. Some studies do not supply public
reads or a complete preprocessing recipe; metadata and warnings explain the gap.

Each [study directory](studies) contains `study.yml`, `download.R` and `process.R`.
JSON syntax is used in YAML files so metadata can be read without another R
dependency. Existing scientific internals and recovered original source workflows
remain readable beside each recipe. To modify preprocessing, use a source checkout:

```sh
Rscript studies/vandeputte2017/download.R --dest /path/to/inputs
Rscript studies/vandeputte2017/process.R --input /path/to/inputs --output /path/to/staging/vandeputte2017
```

Recipes may need tidyverse, readxl, DADA2, reference training data or other
study-specific tools. Their requirements remain in the scripts. The installed
package needs only digest and jsonlite. Routine users do not run recipes.

## phyloseq and miaverse

Select a single counts branch and provide explicitly identified sample/feature
tables. These optional conversions retain every count and align tables by names.
They never guess a branch or merge incompatible taxonomy.

```r
counts <- x$vandeputte2017$counts$original
# Set metadata/taxonomy row names to the sample/feature identifiers first.
ps <- mutt_as_phyloseq(counts)
tse <- mutt_as_miaverse(counts)  # TreeSummarizedExperiment, usable with mia
# Optional: metadata = sample_table, taxonomy = feature_table
```

Install phyloseq or TreeSummarizedExperiment with BiocManager only if needed.
Include microbial-load columns in the supplied sample metadata to carry them into
the converted object.

## Contribute a study

Create `studies/STUDY_ID/{study.yml,download.R,process.R}`. Record DOI/citation,
sequencing/load modalities, accession numbers, exact inputs, availability limits
and processing decisions. Use the standard CLI above; keep scientific decisions
study-specific. Produce a staged `study.rds` and `provenance.json`, validate it,
and build a new immutable release. See [maintainer instructions](docs/RELEASING.md)
and [architecture rules](AGENTS.md).

MUTT provides curated taxonomic data. Inferred functions/pathways and general
downstream microbiome analyses are outside its scope.
