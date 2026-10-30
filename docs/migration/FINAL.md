# Reviewing the MUTT migration

This migration changes how MUTT delivers its studies. The previous package read
study files from a repository checkout and ran their parsers when `mutt()` was
called. The new package retrieves a saved study object from a versioned data
release. Separating the package from its data keeps installation small and lets a
user download just the studies needed for an analysis.

The scientific reference is the output of the previous API at commit `01c648b`.
The migration preserves that output's counts, proportions, taxonomy, load
measurements, metadata, demographics, identifiers, missing values and nested
scientific representations. The only deliberately removed field is the empty
inferred-function field. Existing processed inputs were read; scientific
preprocessing was not rerun. The detailed comparisons and one bounded numerical
precision exception are described in [the scientific audit](AUDIT.md).

## What changes for a user

Users now select study IDs from `mutt_studies()`, then request one or more studies:

```r
mutt_studies()
x <- mutt("vandeputte2017")
y <- mutt(c("vandeputte2017", "contijoch2019"))
```

`mutt()` returns a named list, including for a single study. It verifies the
download's size and SHA-256 checksum, then caches the verified data for reuse.
The installed package needs only `digest` and `jsonlite`; ordinary loading does
not require the old parsers or their preprocessing dependencies.

The default study object preserves the old call with `rawdata=FALSE`,
`align_samples=FALSE` and `functional=FALSE`. A separately stored representation
preserves the old sample-aligned output:

```r
aligned <- mutt("vandeputte2017", representation = "aligned")
```

The default keeps the original sample sets, so counts, metadata and load tables
can contain different samples. The aligned form preserves the previous package's
study-specific alignment rules. These forms retain the study's nested assay,
classifier and cohort branches.
Loading a study does not choose or flatten those branches. Source information,
original input downloads and recipe locations are available separately through
`mutt_source()`, `mutt_download_raw()` and `mutt_recipe()`. Optional phyloseq and
miaverse conversion requires the user to select a counts branch explicitly.
See [the README](../../README.md) for the user workflow.

## Which studies are included

The catalog accounts for 63 source records:

| Records | Migration result |
| --- | --- |
| 40 registered studies | Usable old-API output in both default and aligned forms; both forms are packaged. |
| 17 registered studies | Existing parsers failed in both forms; their errors remain documented and their processed data are unavailable. |
| 6 unregistered source directories | Sources are recorded without inventing a processed output. |

The migration does not repair failed parsers by substituting an older saved
study. Original study inputs remain referenced by their historical archive and
checksums, and the original Git history is retained. Functional-inference code,
resources and tooling have been removed; original authors' source-method notes
remain with the studies.

## What has been verified

Complete-object comparisons cover all 40 usable studies in both forms. The
production loader also reproduced all 80 complete default/aligned objects
exactly, including download verification, extraction and cache reuse, with one
ZIP transfer per study. The committed comparison reports establish the original
capture checks; release validation repeated the checks against the packaged
objects. These are comparisons of complete objects, not just dimensions or a
sample of counts.

The 40 ZIP assets, totaling approximately 286 MiB, have been uploaded to the
[draft data release](https://github.com/Silverman-Lab/mutt/releases),
downloaded again and verified against staging. The draft requires authorized
repository access. Each asset includes provenance describing its inputs,
capture parameters and object checksums.

Public downloads are not live yet. The installed manifest still marks these
assets `staged`. Publication and verification of the public pinned URLs must
precede an `available` manifest; [the release instructions](../RELEASING.md)
describe those remaining steps.

## Inspect a study during review

The usual `mutt()` download examples apply after a study becomes `available`.
For this review, you can inspect the saved data directly: download one familiar
study ZIP from the draft release and read its `study.rds` file. For example,
after downloading `vandeputte2017.zip` into your working directory:

```r
review_dir <- tempfile("mutt-review-")
dir.create(review_dir)
unzip("vandeputte2017.zip", exdir = review_dir)
study <- readRDS(file.path(review_dir, "study.rds"))
names(study)
```

The same ZIP contains `aligned.rds` and `provenance.json`. The first holds the
aligned form; the second records the inputs and checks used for that study.
Inspecting these saved objects does not run the study's preprocessing.

## Where to focus the review

Start with the user examples in [README.md](../../README.md), then follow the
scientific and operational responsibilities separately:

- [AUDIT.md](AUDIT.md) defines the scientific baseline, records study accounting
  and explains the Dreier numerical precision exception.
- [R/mutt.R](../../R/mutt.R), [R/manifest.R](../../R/manifest.R) and
  [R/cache.R](../../R/cache.R) implement study selection, pinned downloads,
  verification and caching. [The fixture tests](../../tests/testthat) exercise
  these paths without downloading real studies.
- [R/source.R](../../R/source.R) and [the study directories](../../studies)
  separate source discovery and maintainer recipes from normal data loading.
  [R/converters.R](../../R/converters.R) requires explicit sample/feature
  alignment for optional conversions.
- [The release instructions](../RELEASING.md) describe staging, archive checks,
  draft verification and publication. Permanent release scripts read staged
  objects and provenance; they do not rerun legacy parsers.

The revised package passed `R CMD check` on R 4.6.1 with Status OK: no errors,
warnings or notes. The installed-package tests passed 953 assertions; the
previous maintainer release-tool run passed five assertions.
Repository-index access was limited during the check. Optional converter object
construction was not fully exercised, and TreeSummarizedExperiment was absent;
explicit identifier/count alignment was tested. The scientific audit describes
what these preservation checks establish and where their coverage ends.
