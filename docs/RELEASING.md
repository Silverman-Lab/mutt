# Maintainer releases

Normal users install only the small R package. No data live in the active Git
source tree. Keep maintainer inputs and staging in gitignored `.build/`.

## One-time v1 capture

The pre-migration baseline is commit `01c648b1f736b1af019d1c0f0231b73de2c21e49`.
Keep its installed package in an isolated R library and materialize its original
non-functional study inputs separately. Never rerun scientific profiles to test
distribution. The local `.build/legacy/` snapshot preserves those files; the
original commit also remains accessible in unchanged Git history.

```sh
TZ=UTC Rscript scripts/migration/capture-current.R OLD_R_LIBRARY LEGACY_STUDIES PARSER_INDEX .build/capture
Rscript scripts/migration/stage-current.R .build/capture/default .build/capture/aligned .build/data-v1.0.0
```

Capture invokes the old API on stored taxonomic tables and explicitly blocks
DADA2 classification/denoising. Errors are recorded rather than repaired with new
science or substituted historical objects. Only the deliberately removed empty
functional-inference field is stripped. Compare serialization and scientific
objects exactly. The committed comparison report and per-study provenance record
the result. A historical 37-study RData was investigated, but the user's clarified
requirement selects current API outputs instead. `stage-v1.R` is retained solely
as the historical extraction record; do not use it to select v1 release data.

After initial migration these scripts are historical, unsupported ongoing tools.
New studies enter staging through their normal study recipes.

Some old parsers require the `dada2` namespace even to read stored tables.
Install that parsing dependency if needed, while retaining the capture script's
classification and denoising guards. See `docs/migration/FINAL.md` for the
approved, checksum-verified Dreier scale recovery used on the laptop.

## Permanent staged contract

Each staged directory is named by study ID and contains:

- `study.rds`: the complete default scientific object;
- `provenance.json`: study ID, schema_version, data_version, selection notes,
  metadata and `study_rds_sha256`;
- optional named RDS forms declared in `variants` with file and SHA-256, e.g.
  `aligned.rds` for the previously sample-aligned output.

The default object has counts, proportions, tax, scale and metadata fields. Nested
assays and other legitimate fields are allowed. Do not reconstruct an RDS for
cosmetic reasons. Validate loaded scientific contents without normalizing values.

```sh
Rscript scripts/validate-study.R .build/data-v1.0.0/vandeputte2017
Rscript scripts/package-study.R .build/data-v1.0.0/vandeputte2017 .build/one-asset data-v1.0.0
Rscript scripts/validate-asset.R .build/one-asset/vandeputte2017.zip .build/data-v1.0.0/vandeputte2017
Rscript scripts/build-release.R --input .build/data-v1.0.0 --tag data-v1.0.0 --output .build/release-v1
```

These permanent tools read only staging/provenance, with no legacy-path logic.
They reject overwrites, check object/provenance/variant checksums, create one ZIP
per study, verify archive layout, load objects, compare scientific contents and
compare extracted bytes to staging. ZIP packaging requires a `zip` command.
The candidate manifest is generated in the asset directory with status `staged`.

## Draft upload and download verification

Authenticate `gh` with release write access to Silverman-Lab/mutt and push the
reviewed source commits. An authenticated successful release listing is required
before any mutation. A published tag is never modified. Resumed draft uploads do
not clobber existing assets; those must verify exactly against staging.

```sh
Rscript scripts/publish-draft.R .build/data-v1.0.0 .build/release-v1
```

This creates a draft release, uploads missing assets, downloads **every** asset,
checks size/SHA-256/archive integrity, loads all study forms and compares them to
the staged source. `draft-verification.json` is written only after all checks pass.
A failed authentication or transfer is not a successful release. The package
must never claim staged assets are already public.

Publication is a separate explicit maintainer action. After publication, run
`verify-public.R` to verify the unauthenticated pinned URLs and generate an
available manifest. Merge those generated entries into the package manifest,
retaining explicit unavailable records. Never mark an asset available before
that public verification. Never use latest or overwrite a published release;
corrections use a new data tag and an updated pinned package manifest.

## Validation

Run local fixture tests and build/check the source package. CI must not run
migration, download study data/FASTQs or build a data release. Optional Bioconductor
converters are not dependencies of the ordinary download API.

The maintainer-only release fixture test lives outside the installed package:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE); source("tests/testthat/helper-release.R"); testthat::test_file("scripts/tests/test-release-tools.R", stop_on_failure=TRUE)'
```

It packages only tiny local fixtures. Installed-package tests must not depend on
maintainer scripts excluded from the source package.
