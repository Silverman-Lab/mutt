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
remain in a gitignored recovery snapshot; complete historical study objects are
not substituted for current data. The bounded Dreier scale recovery is documented
below.

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

## Local Dreier recovery

A laptop recapture of Dreier differed from the cloud reference in 150 values:
six load-summary columns across 25 samples. Counts, proportions, taxonomy,
metadata, demographics, identifiers and missingness were identical. The maximum
relative difference was `2.0302307561987318e-14`; maximum absolute differences in
`mean_copies` and `sd_copies` were `9.1968104243278503e-09` and
`8.4401108324527740e-09`, respectively. Local R 4.6.1 versus cloud R 4.5.3 is a
plausible source of numerical variation, but the cause is not fully established.

The recovery uses only the historical `scale` field from Git blob
`5ef175e49c7eec8b4725cf33ac5aa2bf087bbbc8` with the otherwise unchanged local
current-API object. The complete reconstructed default object exactly matches
the recorded cloud scientific SHA-256:
`4951997ad8467ab87320b35a69e38db2a4792ef2dff3206bd4094be17a45effc`.
The full historical study object does **not** match and is not substituted.
This preserves the exact baseline; it does not relax equality to a numerical
tolerance. The local audit is retained under
`.build/laptop-release/dreier-audit/` in `reference-verification.json`,
`component-comparison.csv` and `scale-difference-summary.csv`. Audit verification
alone does not establish release staging or draft upload/download verification.

## Operational state

40 independently downloadable ZIPs have been built locally for `data-v1.0.0`.
The permanent builder validates canonical staging, object/provenance/variant
checksums, archive layout, scientific object loading and equality with staging.
It does not know legacy paths. The draft uploader never clobbers assets and
downloads every uploaded asset for verification. Public URL verification is a
separate explicit step after publication.

The cloud environment's gh request returned Forbidden, preventing upload there.
The laptop continuation on 2026-10-29 captured 40 successful studies and the same
17 parser failures in each form. Pereira loaded after installing its existing
`dada2` parsing dependency. With UTC capture and the bounded Dreier recovery,
all 40 default scientific checksums match the committed cloud baseline. All 40
local ZIPs passed the permanent builder; the production loader reproduced all
80 complete default/aligned objects exactly using one ZIP transfer per study.
Each asset records the local capture environment; Dreier also records the raw
local hashes, recovered baseline hashes and the approved recovery scope.

Local evidence lives in `.build/laptop-release/capture-final/`,
`loader-verification/` and `release-v1/`. Remote upload/download verification is
recorded separately in `release-v1/draft-verification.json` only after success.
The installed manifest remains staged/unavailable until publication and public
URL verification. Do not mark staged assets available without those checks.
See `docs/RELEASING.md` for the concrete commands.

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
Final package-check and fixture logs accompany the deliverables. The local
R 4.6.1 source build/check completed with zero errors, zero check warnings and
one NOTE about unrecognized NEWS entries; 953 installed fixture assertions and
five maintainer release-tool assertions passed. Repository-index access was
unavailable during the check, and TreeSummarizedExperiment was not installed.

A real study-recipe smoke check also reconstructs Vandeputte 2021 from stored
processed inputs and compares the entire aligned object exactly against the old
API. JSON null demographic fields are preserved as NULL, not empty lists.
