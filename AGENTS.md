# MUTT architecture

MUTT is a curated, versioned R library for harmonized microbiome datasets,
particularly taxonomic profiles paired with microbial-load measurements.
Users install the small `mutt` R package and request only needed studies with
`mutt_studies()` and `mutt(c("vandeputte2017", "contijoch2019"))`. They should
never need a full data clone, LFS, Docker, DADA2 or MetaPhlAn for ordinary use.

## Storage and reproducibility

Git stores code, study metadata, download/processing scripts, small fixture tests,
documentation and pinned manifests. GitHub Releases store independently
downloadable, immutable processed study ZIPs. Raw reads remain in authoritative
repositories such as SRA/ENA. The R user cache contains only requested processed
datasets and explicitly requested source inputs. `.build/` is gitignored and
contains maintainer staging, assets and one-time baseline snapshots.

Each asset contains `study.rds`, `provenance.json`, and optionally named RDS
representations declared with checksums in provenance. Scientific objects have
counts, proportions, tax, scale and metadata, possibly nested by assay, classifier,
cohort or specimen type. Preserve legitimate branches, missing values, units and
identifiers; never normalize scientific values during distribution. `aligned.rds`
preserves the old API's sample-aligned form. Default uses the old unaligned form.

The installed manifest is pinned to explicit `data-vX.Y.Z` tags, sizes and SHA-256.
Never resolve latest. Never overwrite published assets: corrections need a new
release. Staged/draft assets must be listed honestly and cannot be treated as
publicly downloadable. Cache extraction is verified in staging and published
with a same-filesystem rename; failed refresh leaves the previous copy usable.

## Study contract

Every study directory has `study.yml`, `download.R`, `process.R`:

    Rscript studies/STUDY/download.R --dest /inputs [--reads]
    Rscript studies/STUDY/process.R --input /inputs --output /staging/STUDY

`study.yml` uses JSON syntax (valid YAML 1.2). Sources vary: some studies provide
only processed tables, some have raw reads, some lack a complete recipe. Say so
explicitly. Missing public data generate useful warnings; never invent URLs,
accessions, checksums or complete processing recipes. Historical source archives
are immutable fallbacks; prefer authoritative sources where exact mappings exist.

Standardize the external script interface; preserve study-specific internals.
Shared helpers may handle simple retrieval, directories, checksums and validation.
Do not centralize heterogeneous scientific processing into one giant pipeline.
Recipe dependencies are maintainer-only and do not belong in package Imports.

## V1 baseline and one-time migration

The pre-migration manuscript-era artifacts are the v1 scientific baseline.
The user's clarified requirement is exact equality with **current mutt at
01c648b**, including counts, metadata, taxon names and all branches, allowing
interface edits. Do not unnecessarily regenerate scientific profiles.
`scripts/migration/` isolates old paths and parser knowledge. It captures current
API outputs from stored processed tables without denoising/profiling, strips only
the deliberately removed empty functional field, records original artifact hashes,
and copies RDS bytes into clean staging. Old broken parsers remain explicitly
unavailable; do not silently substitute an earlier historical object.

The separate 37-study historical RData was investigated and retained locally for
recovery; it is not the selected baseline. Permanent tools in `scripts/release.R`
know only staged objects and provenance, never old paths. They validate/package,
generate candidate manifests, create draft releases with gh, upload without
clobbering, download and compare every asset to staging. Only after public
publication and independent URL verification may manifest statuses be available.

## Scope and maintenance

PICRUSt2, FAPROTAX, inferred functional profiles/pathways and downstream biological
analysis are intentionally outside MUTT. DADA2, RDP and MetaPhlAn may construct
study taxonomic data on the maintainer side. Do not retain functional inference
containers, HPC jobs or bundles. Routine CI checks only package code and tiny
local fixtures, cancels superseded runs, and never rebuilds the database.

Do not introduce speculative subsystems. Before major architecture changes,
explain the scientific requirement, why it belongs here, why simpler code is
insufficient, and maintenance, storage and reproducibility implications.
Never rewrite Git history as part of distribution cleanup.
