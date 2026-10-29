# V1 migration result

The scientific reference is current mutt at
`01c648b1f736b1af019d1c0f0231b73de2c21e49`, following the user's clarification.
Only its empty inferred-function field is removed. Counts, proportions, metadata,
taxon names, taxonomic branches, load values, sample identifiers, other scientific
fields and demographics are preserved without denoising or profiling.

## Data accounting

Both default and aligned forms were captured for all 57 registered studies.
40 succeeded in each form; 17 existing parsers failed in each form. Full warnings
and errors are in `current-default-parser-report.csv` and
`current-aligned-parser-report.csv`. Six unregistered legacy directories remain
source records without invented processed outputs. The catalog consequently has
63 source records and 40 staged processed studies. Historical 37-study objects
remain in a gitignored recovery snapshot and are not substituted for current data.

`current-object-comparison.csv` records byte-for-byte and complete object equality
between staged objects and captured old API results. `loader-comparison.csv`
records 40/40 default and 40/40 aligned exact end-to-end comparisons and equality through the production downloader, verification,
extraction and cache with local ZIP transfers replacing network access. Each study
is retrieved once for its two stored forms. No comparison is limited to dimensions
or a sample of counts. Every asset contains original input hashes and capture
parameters in its provenance.

The active tree contains code, metadata, recipes, tests, documentation and pinned
manifests. It contains no study data or LFS attributes. The original sources remain
in `.build/legacy`, immutable archive references, and unchanged Git history. The
installed R package imports only digest and jsonlite. Recipe dependencies remain
maintainer-only. Normal users request only selected studies.

## Operational state

40 independently downloadable ZIPs have been built locally for `data-v1.0.0`.
The permanent builder validates canonical staging, object/provenance/variant
checksums, archive layout, scientific object loading and equality with staging.
It does not know legacy paths. The draft uploader never clobbers assets and
downloads every uploaded asset for verification. Public URL verification is a
separate explicit step after publication.

The environment's gh request returned Forbidden; no draft release or asset upload
has succeeded. The manifest says staged for local assets and unavailable for
other records. Ordinary first downloads cannot succeed until upload, publication
and public URL verification are completed. Do not mark staged assets available
without those checks. See `docs/RELEASING.md` for the concrete commands.

## Scope and validation

Functional inference implementation, resources, exports, tests, HPC products and
containers were removed. Remaining search hits document that decision or quote
original author methods. In particular, Wagner's source notes describe original
authors' 16S copy-number correction; that historical context and its counts were
preserved. A recovered Regalado script contains non-R console/shell material;
it is retained verbatim as `correct_16S_load.R.txt`, not presented as executable.

Fixture tests cover selected/multiple downloads, checksum/size failures, partial
cleanup, cache reuse/corruption/refresh/removal, aligned forms, manifest pinning,
archive members, staged release tooling and explicit conversion-table alignment.
Routine CI checks small fixtures across three operating systems and cancels
superseded runs; it never rebuilds studies or downloads FASTQs. Actual phyloseq
and TreeSummarizedExperiment construction has not been exercised here because
those optional packages are not installed; identifier/count alignment is tested.

Local data verification does not replace remote upload/download verification.
Final package-check and fixture logs accompany the deliverables.

A real study-recipe smoke check also reconstructs Vandeputte 2021 from stored
processed inputs and compares the entire aligned object exactly against the old
API. JSON null demographic fields are preserved as NULL, not empty lists.
