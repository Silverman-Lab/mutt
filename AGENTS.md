# Working on MUTT

MUTT is a curated R library of microbiome taxonomic profiles paired with
microbial-load measurements. The R package now retrieves processed study objects
from separately versioned data releases. Users list studies with `mutt_studies()`
and request only those they need with `mutt(c("vandeputte2017", "contijoch2019"))`.
Ordinary use must not require a full data clone, LFS, Docker, DADA2 or MetaPhlAn.

## Code, data and versions

Git stores package code, study metadata, download and processing recipes, small
fixture tests, documentation and pinned manifests. GitHub Releases store one ZIP
per processed study. Raw reads remain in authoritative repositories such as
SRA/ENA. Keep maintainer staging, assets and historical recovery snapshots in
gitignored `.build/`; do not add study data to the active source tree.

The installed manifest is the catalog that connects study IDs to release assets.
It pins explicit `data-vX.Y.Z` tags, URLs, sizes and SHA-256 checksums. Never
resolve `latest` or overwrite published assets. Corrections require a new data
tag and an updated manifest. A staged or draft asset must not be advertised as
publicly downloadable. Follow [docs/RELEASING.md](docs/RELEASING.md) for the
current release status and maintainer workflow.

## Preserve scientific objects

Each study ZIP contains `study.rds`, `provenance.json` and any named RDS forms
declared with checksums in provenance. Scientific objects contain counts,
proportions, tax, scale and metadata. These fields may be nested by assay,
classifier, cohort or specimen type. Preserve all legitimate branches, missing
values, units, sample identifiers and other scientific fields. Distribution
must not normalize or otherwise change scientific values.

The default form preserves the old API's `rawdata=FALSE, align_samples=FALSE`
output. The optional `aligned.rds` preserves its `align_samples=TRUE` output.
Do not align or select an assay implicitly when returning the default object.

For v1, the scientific baseline is the API output at commit
`01c648b1f736b1af019d1c0f0231b73de2c21e49`, with only the deliberately removed
empty `function` field omitted. Require exact equality of complete scientific
objects. The separate historical 37-study RData is a recovery reference and must
not replace that baseline. A bounded recovery of Dreier's `scale` field preserved
exact baseline equality; its scope and checks are recorded in
[docs/migration/FINAL.md](docs/migration/FINAL.md).

The one-time scripts in `scripts/migration/` isolate old paths and parser
knowledge. They read stored processed tables with classification and denoising
blocked, record input hashes and capture parameters, and copy captured RDS bytes
into staging. Use UTC for historical capture. Existing broken parsers remain
explicitly unavailable; do not invent outputs or silently substitute historical
study objects. Future studies use their normal recipes and the permanent staged
release tools, which must not acquire legacy-path logic.

## Downloads, validation and cache safety

The R user cache contains only requested processed studies and explicitly
requested source inputs. Verify asset sizes and checksums, expected archive
members, provenance and every stored form before accepting a study. Reject
unexpected archive paths and undeclared files. Keep extraction in a temporary
directory and publish the verified cache with a same-filesystem rename. A failed
refresh must leave the previous cached copy usable.

Permanent tools in `scripts/release.R` read staged objects and provenance. They
validate and package studies, generate candidate manifests, upload drafts without
clobbering assets, and download every asset to compare it with staging. Only
publication followed by independent public URL verification permits status
`available` in the package manifest.

## Study recipes

Each study directory contains `study.yml`, `download.R` and `process.R`, with
this common external interface:

```sh
Rscript studies/STUDY/download.R --dest /inputs [--reads]
Rscript studies/STUDY/process.R --input /inputs --output /staging/STUDY
```

`study.yml` uses JSON syntax, which is valid YAML 1.2. Source coverage varies:
some studies provide processed tables, some provide raw reads, and some lack a
complete recipe. Record these differences explicitly. Warn clearly about missing
public inputs; never invent URLs, accessions, checksums or complete recipes.
Historical source archives are immutable fallbacks. Prefer authoritative sources
where exact mappings exist.

Keep the common interface while preserving study-specific scientific processing.
Shared helpers may handle retrieval, directories, checksums and validation; do
not combine heterogeneous processing into one large pipeline. Recipe dependencies
are maintainer-only and do not belong in the installed package's Imports.

## Scope and review

PICRUSt2, FAPROTAX, inferred functional profiles and pathways, and downstream
biological analysis are outside MUTT's scope. Do not restore their implementation,
containers, HPC jobs or bundles. DADA2, RDP and MetaPhlAn may construct taxonomic
data in maintainer recipes. Preserve original scientific source counts and notes
about author methods, including historical copy-number correction.

Routine CI checks package code and tiny local fixtures, cancels superseded runs,
and never rebuilds the study database or downloads study data or FASTQs. Keep
installed-package tests independent of maintainer scripts excluded from the
source package.

Before a major architecture change, explain the scientific requirement, why it
belongs in MUTT, why simpler code is insufficient, and the maintenance, storage
and reproducibility implications. Never rewrite Git history as part of
distribution cleanup.
