# mutt

MUTT brings together microbiome studies that pair taxonomic profiles with
microbial-load measurements. It provides processed study data, source information
and the study-specific scripts used to assemble those data.

This update makes the package easier to install and use. The R package contains
the study catalog and loading functions. Processed data are stored separately,
with one downloadable ZIP per study. You request the studies you need, and MUTT
keeps a checked copy on your computer for later use.

The existing scientific data are preserved. Changing how a study is downloaded
does not rerun its sequencing analysis or choose among its scientific results.

## Status of this update

The 40 study ZIPs have been uploaded to the [draft data release](https://github.com/Silverman-Lab/mutt/releases)
and downloaded again for verification. All 80 default and aligned study objects
passed complete-object comparisons. The release remains a draft for review, so
public study downloads are not enabled in the package yet.

For an overview of the changes and suggested review order, start with the
[review guide](docs/migration/FINAL.md). The [scientific audit](docs/migration/AUDIT.md)
contains the detailed preservation checks.

## Install and explore the catalog

While reviewing this update, install its branch:

```r
install.packages("remotes")
remotes::install_github("Silverman-Lab/mutt", ref = "codex/mutt-v1-release")
library(mutt)

studies <- mutt_studies()
studies[, c("study_id", "display_name", "status")]
```

The catalog uses short study IDs, such as `vandeputte2017`. It also records
whether processed data can be downloaded:

| Status | Meaning |
| --- | --- |
| `available` | The published download has been verified and can be requested. |
| `staged` | The processed data are prepared for release; public downloads are not enabled yet. |
| `unavailable` | This release has no usable processed object for the study. Source information may still be available. |

Listing studies and inspecting their source information does not download study
data. The current catalog contains 40 prepared studies, 17 studies whose previous
parsers failed, and six additional source records without processed outputs.

## Load a processed study

Once a study is marked `available`, request it by ID:

```r
x <- mutt("vandeputte2017")
study <- x$vandeputte2017
names(study)

# Request several studies together.
y <- mutt(c("vandeputte2017", "contijoch2019", "morton2019"))
```

`mutt()` returns a named list, including when you request a single study. The
study objects retain their existing fields:

| Field | Contents |
| --- | --- |
| `counts` | Sequence count tables. |
| `proportions` | Relative abundance tables. |
| `tax` | Taxonomic annotations. |
| `scale` | Microbial-load measurements, with their original names and units. |
| `metadata` | Sample information. |

Additional fields, including demographics, are retained where present. A field
can contain several named tables for different assays, classifiers or cohorts.
Those alternatives remain separate so you can select the appropriate one for
your analysis.

### Sample alignment

The default form retains the sample sets in the previous package's unaligned
output. Counts, metadata and load tables can therefore contain different samples.
To retrieve the previous package's sample-aligned form, use:

```r
aligned <- mutt("vandeputte2017", representation = "aligned")
```

This preserves each study's existing alignment behavior. Check sample identifiers
when combining its tables. The [scientific audit](docs/migration/AUDIT.md)
explains how both forms were compared with the previous package.

### Reuse a downloaded study

Each download is checked against its recorded size and SHA-256 checksum, a
fingerprint of the file's exact contents. The checked data are kept in R's user
cache, so repeated requests reuse the saved copy.

```r
x <- mutt("vandeputte2017")                 # reuse the checked copy
x <- mutt("vandeputte2017", refresh = TRUE) # download and check again
mutt_cache_remove("vandeputte2017")        # remove the saved copy
```

Routine loading needs neither a repository checkout nor scientific preprocessing
software. The package's required R dependencies are `digest` and `jsonlite`.

## Inspect sources and processing scripts

```r
mutt_source("vandeputte2017") # source notes, accessions and input information
mutt_recipe("vandeputte2017") # links to the study's processing scripts
```

If you need the original inputs, request them separately:

```r
mutt_download_raw("vandeputte2017", dest = "inputs")
# Also request public sequencing reads, when available:
# mutt_download_raw("vandeputte2017", dest = "inputs", reads = TRUE)
```

Despite its name, `mutt_download_raw()` can retrieve published tables and stored
processed inputs as well as raw sequencing reads. Reads are fetched from ENA only
when `reads = TRUE`. Some studies lack public reads or complete processing
instructions; their source information describes these limits.

Each [study directory](studies) contains metadata in `study.yml`, a `download.R`
script and a `process.R` script. These scripts are for maintaining or extending
the collection. Their study-specific dependencies, such as DADA2 or readxl, are
separate from the dependencies needed to load released data.

## Use phyloseq or miaverse

Choose a single counts table before converting it:

```r
counts <- study$counts$original
ps <- mutt_as_phyloseq(counts)
tse <- mutt_as_miaverse(counts)
```

Counts must have samples in rows and features in columns. You can also supply
`metadata` and `taxonomy` tables with the corresponding identifiers as row names.
Include microbial-load columns in the sample metadata if you want them in the
converted object. Conversion keeps the selected counts and matches the other
tables by identifier.

Install the optional `phyloseq` or `TreeSummarizedExperiment` package with
BiocManager when using its converter.

## Review or contribute

The [review guide](docs/migration/FINAL.md) introduces the package changes and
points to the relevant code. The [scientific audit](docs/migration/AUDIT.md)
records data coverage, exact comparisons and the Dreier numerical precision
note. [Maintainer instructions](docs/RELEASING.md) describe adding a study and
preparing a data release.

MUTT focuses on curated taxonomic profiles and microbial-load measurements.
Functional inference and downstream biological analyses are outside this scope.
