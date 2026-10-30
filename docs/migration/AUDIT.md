# Scientific migration audit

This audit explains what “preserved” means in the MUTT migration and identifies
the evidence for that claim. The migration stores the previous API's study
objects in downloadable assets. It changes data delivery while retaining the
scientific contents of those objects.

## Reference and preservation rule

The reference is package `0.1.0.9000` at commit
`01c648b1f736b1af019d1c0f0231b73de2c21e49`. The previous public call was:

```r
mutt(studies = NULL, base_directory = NULL, rawdata = FALSE,
     align_samples = FALSE, save_to = NULL, verbose = FALSE,
     functional = FALSE)
```

The default reference uses `rawdata=FALSE`, `align_samples=FALSE` and
`functional=FALSE`. A second reference uses `align_samples=TRUE`. Captures read
the stored processed inputs through the previous API. Only its deliberately
removed empty `function` field is stripped. Counts, proportions, taxon names,
taxonomy, load measurements, sample identifiers, metadata, demographics, missing
values, attributes and nested scientific representations otherwise remain part
of the complete-object comparison.

The reference is the current API at that commit, rather than a reconstruction
of a manuscript dataset. Manuscript and supplementary descriptions provide
context for original/reprocessed, RDP16, MetaPhlAn and mOTU-derived branches;
they do not establish that one branch should replace the others. Legitimate
alternate representations are retained.

No denoising, taxonomy classification, MetaPhlAn profiling or other scientific
preprocessing was rerun. The one-time capture script guards DADA2 classification
and denoising entry points so that a missing profile fails explicitly. Some
parsers still require the DADA2 namespace to read existing tables; that
dependency does not establish that scientific preprocessing ran. Existing
parser failures were recorded rather than repaired with new profiles or older
complete study objects.

## Source preservation and study accounting

[baseline-files.txt](baseline-files.txt) records the tracked inventory.
[baseline-studies.json](baseline-studies.json) records registered parsers,
artifacts, source notes, URLs, accessions and study modalities.
[baseline-artifacts.csv](baseline-artifacts.csv) contains hashes, sizes and LFS
status for 847 original study artifacts. All 847 original artifact hashes match
the preserved source inputs. Historical archive references and unchanged Git
history preserve access to those inputs outside the small package source tree.
An accession reference alone does not establish that every original input is
publicly obtainable; availability limits remain in the source metadata.

The old parser index registered 57 studies. Both representations were attempted
for every registered study: 40 succeeded and 17 failed in each form. The
[default parser report](current-default-parser-report.csv) and
[aligned parser report](current-aligned-parser-report.csv) retain the full
warnings and errors. The later release capture reproduced the same success and
failure accounting.

Six additional legacy directories were not registered with the old API.
[unregistered-studies.json](unregistered-studies.json) and
[unregistered-source-inventory.json](unregistered-source-inventory.json) record
them as sources without selected harmonized outputs. The catalog therefore has
63 source records: 40 with packaged processed data, 17 with existing parser
failures and six unregistered sources. No processed output was invented for the
last two groups.

## Complete-object and loader comparisons

[current-object-comparison.csv](current-object-comparison.csv) records exact
comparison results for both forms of all 40 usable studies, along with
serialized-file and scientific-object SHA-256 values for the default form. These committed results
describe the original reference capture. They establish agreement of complete
objects rather than structural agreement alone.

[loader-comparison.csv](loader-comparison.csv) records exact comparisons for all
40 default and all 40 aligned objects through the production downloader,
checksum verification, extraction and cache paths. Local ZIP transfers supplied
the archives for this test. Each study required one transfer for its two stored
forms. The release validation repeated complete-object loader comparisons for
all 80 forms and confirmed that all 40 default scientific hashes matched the
committed reference, with the bounded Dreier reconstruction described below.

Every release asset contains `provenance.json` beside the study RDS files. It
records original input hashes, capture parameters, environment information and
checksums for the default and aligned forms. Dreier's provenance additionally
distinguishes the raw recapture hashes from the reconstructed reference hashes
and records the scope of the scale-field recovery. File checksums verify the
packaged bytes; complete-object comparisons verify their scientific contents.

## Dreier: cross-environment numerical precision

The Dreier recapture differed from the reference in 150 values: six
load-summary columns across 25 samples. Counts, proportions, taxonomy, metadata,
demographics, identifiers and missingness were identical. The observed
differences were small but nonzero:

| Comparison | Maximum difference |
| --- | --- |
| Relative difference | `2.0302307561987318e-14` |
| Absolute difference in `mean_copies` | `9.1968104243278503e-09` |
| Absolute difference in `sd_copies` | `8.4401108324527740e-09` |

The reference and recapture used R 4.5.3 and R 4.6.1, respectively. Numerical
variation between environments is a plausible explanation; the cause was not
fully established. The unadjusted recapture therefore differs from the reference at this precision.

The reconstruction replaces only the `scale` field, using the field recovered
from historical Git blob `5ef175e49c7eec8b4725cf33ac5aa2bf087bbbc8`. All other
fields remain those of the current-API recapture. The complete reconstructed
default object matches the recorded reference scientific SHA-256 exactly:

```text
4951997ad8467ab87320b35a69e38db2a4792ef2dff3206bd4094be17a45effc
```

The recaptured default and aligned scale fields were exactly identical before
recovery, so the aligned form received the same verified scale field. The
release checks then compared the complete aligned object with this recovered
reference. There is no separately committed scientific checksum for the original
aligned reference; the full-object checksum proof above applies to the default
form. Checks of the recovered objects through staging and loading remain exact.

The historical blob contains `datasetsfromrepo.RData`, a saved collection of 37
five-field study objects from commit `53a5406`, deleted in `ebbdc08`. A complete
historical Dreier object does not match the selected current-API reference and
was not substituted. This bounded recovery supplies only the reference scale
field; it does not select the historical collection as the scientific baseline
or repair any failed parser.

## Scientific scope and recipe checks

The migration removes inferred-function implementation, resources, exports,
tests, HPC products and containers. Some original source notes describe methods
that used functional tools. Wagner's authors, for example, describe PICRUSt2
16S copy-number correction. That source-method note and the existing counts
remain intact; they do not introduce a new inference step into this migration.

The recovered Regalado source file contains non-R console/shell material. It is
preserved verbatim as `correct_16S_load.R.txt` for inspection, rather than
presented as executable R. The existing Regalado parser failure remains in the
failure reports.

A Vandeputte 2021 recipe smoke check reconstructed a study from stored processed
inputs and matched the complete old aligned object exactly. JSON null
demographic fields remained R `NULL`, rather than empty lists. This checks a
real recipe without claiming that every study's scientific preprocessing was
rerun.

## Release and package validation

The permanent release builder reads only staged objects and provenance. It
checks object and variant checksums, archive layout, extracted bytes and loaded
scientific contents. It contains no legacy-path selection logic. The draft
uploader verifies existing assets without clobbering them, downloads every
uploaded asset and compares it with staging.

All 40 study ZIPs, approximately 286 MiB, have passed those checks after upload
and download from the
[draft release](https://github.com/Silverman-Lab/mutt/releases).
The release tools produce `draft-verification.json` only after successful remote
verification; this maintainer report is distinct from the committed original
capture and loader CSVs. Public URL verification remains a separate step after
publication. The installed manifest therefore still marks the assets `staged`.

The revised source package passed `R CMD check` on R 4.6.1 with Status OK:
zero errors, warnings or notes. The installed-package tests passed 953 assertions.
The separate maintainer release-tool run passed five assertions. Fixtures cover
selected and multiple downloads, checksum/size failures, partial cleanup,
cache reuse/corruption/refresh/removal, stored aligned forms, manifest pinning,
archive members, release tools and explicit conversion-table alignment. Routine
CI uses small fixtures and does not rebuild studies or download FASTQs.

Repository-index access was limited during the package check.
TreeSummarizedExperiment was not installed, and actual optional converter object
construction was not fully exercised; identifier/count alignment was tested.
These checks establish the tested implementation and its limits; they do not
independently revalidate the original studies' scientific methods.
