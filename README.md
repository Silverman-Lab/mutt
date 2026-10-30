# mutt

MUTT provides processed microbiome study data that pair taxonomic profiles with
microbial-load measurements. You choose the studies you need, and MUTT downloads
them and saves a copy on your computer for later use.

You can load existing study data directly from R. You do not need to download the
whole project or run read-processing tools.

## Install and find studies

```r
install.packages("remotes")
remotes::install_github("Silverman-Lab/mutt")
library(mutt)

studies <- mutt_studies()
studies[, c("study_id", "display_name", "status")]
```

Use the short study IDs, such as `vandeputte2017`, to request data. Only studies
marked `available` can be downloaded.

| Status | Meaning |
| --- | --- |
| `available` | Processed study data can be downloaded. |
| `staged` | Processed study data are prepared but cannot yet be downloaded. |
| `unavailable` | No processed data are provided for this study in this release. Source information may still be available. |

Listing studies and inspecting their source information does not download study
data. The loading examples below require an `available` study. If none are
available, you can still use the source-information functions.

## Load study data

```r
x <- mutt("vandeputte2017")
study <- x$vandeputte2017
names(study)

# Request several studies together.
y <- mutt(c("vandeputte2017", "contijoch2019", "morton2019"))
```

`mutt()` returns a named list, even when you request one study. Each study contains
the following fields:

| Field | Contents |
| --- | --- |
| `counts` | Sequence count tables. |
| `proportions` | Relative abundance tables. |
| `tax` | Taxonomic annotations. |
| `scale` | Microbial-load measurements, with their original names and units. |
| `metadata` | Sample information. |

Additional fields, including demographics, are included where present. A field
can contain several named tables for different assays, classifiers or cohorts.
For example, `names(study$counts)` shows the count tables you can choose from.
Select the table you need for your analysis.

### Which samples are included?

By default, MUTT loads the saved study tables with their original sample groups.
The counts, sample information and microbial-load tables may contain different
samples.

To work with counts and load measurements from the same samples, some studies
also provide a saved version prepared by matching sample IDs. Request it with:

```r
aligned <- mutt("vandeputte2017", representation = "aligned")
```

Here, “aligned” refers to matching sample IDs between tables, not aligning DNA
sequences. For example, if counts contain samples A, B and C but load measurements
exist only for A and C, this version may limit the counts to A and C. The rules
vary by study: they do not guarantee that every table contains the same samples
in the same order. Check sample IDs when combining tables, rather than matching
rows by position.

### Reuse a downloaded study

MUTT checks each download and keeps a local copy, called a cache. Later requests
reuse that copy.

```r
x <- mutt("vandeputte2017")                 # reuse the saved copy
x <- mutt("vandeputte2017", refresh = TRUE) # download and check again
mutt_cache_remove("vandeputte2017")        # remove the saved copy
```

## Inspect sources and processing scripts

```r
mutt_source("vandeputte2017") # source notes, accessions and input information
mutt_recipe("vandeputte2017") # links to the study's processing scripts
```

If you need the original inputs, download them separately:

```r
mutt_download_raw("vandeputte2017", dest = "inputs")
# Also request public sequencing reads, when available:
# mutt_download_raw("vandeputte2017", dest = "inputs", reads = TRUE)
```

Despite its name, `mutt_download_raw()` can retrieve published tables and stored
processed inputs as well as raw sequencing reads. Reads are downloaded from ENA
only when `reads = TRUE`. Source notes describe missing public inputs and gaps in
processing instructions.

The scripts in each [study directory](studies) are for maintaining or extending
the collection. See the [maintainer instructions](docs/RELEASING.md) for adding
studies and preparing data releases.

## Use phyloseq or miaverse

Install the optional `phyloseq` or `TreeSummarizedExperiment` package with
BiocManager before using its converter. For `vandeputte2017`, the original
counts include results from more than one taxonomic classification. This example
selects the table named `rdp19`:

```r
counts <- study$counts$original$rdp19
ps <- mutt_as_phyloseq(counts)
tse <- mutt_as_miaverse(counts)
```

Counts must have samples in rows and features, such as ASVs or taxa, in columns. You can also supply
`metadata` and `taxonomy` tables with sample IDs and feature IDs, respectively,
as row names. These tables must cover all IDs in the selected counts table; the
converter puts their rows in the corresponding order and keeps all selected
counts. Include microbial-load columns in the sample metadata if you want them
in the converted object.
