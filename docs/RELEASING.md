# Maintaining study data releases

MUTT now distributes its R package and its processed study data separately. The
package contains the loader, study catalog and links to study recipes. Each
processed study is stored in its own ZIP on GitHub Releases, so users download
only the studies they request. Installing or using the package does not run study preprocessing.

The package's **manifest**, `inst/extdata/studies.json`, connects these two parts.
It records each study's availability, an explicit data release tag, its download
URL, file size and SHA-256 checksum. A checksum identifies a file's exact bytes;
the loader checks it before accepting a download. Package versions and data
versions are independent: `data-v1.0.0` names a data release, not an R package
version.

## Current release status

All 40 usable studies have been packaged, uploaded to the `data-v1.0.0` draft,
and downloaded again for verification against their staged sources. This work is
complete. The draft is still unpublished, so the package manifest correctly
keeps these studies marked `staged`. The 17 studies whose existing parsers failed
remain explicitly unavailable; six additional legacy directories have source
records without processed outputs.

The remaining release steps are publication, verification of the public download
URLs, and review of the resulting available manifest. Those steps are described
below. The scientific comparison evidence is in [the migration result](migration/FINAL.md).

## Preparing a future data release

A release begins with **staging**: a local directory containing the exact files
that will go into each study ZIP. Keep staging and release output under the
gitignored `.build/` directory. New studies enter staging through their normal
study recipes; the historical capture scripts are not part of this workflow.

Each staging subdirectory is named by study ID and contains:

- `study.rds`: the complete default scientific object;
- `provenance.json`: the study ID, schema and data versions, selection notes,
  metadata and the checksum of `study.rds`;
- optional named RDS forms declared in the provenance's `variants` entries,
  including a filename and checksum for each form.

The default object preserves the old API's output with `rawdata=FALSE` and
`align_samples=FALSE`. The optional `aligned.rds` preserves its
`align_samples=TRUE` output. Both forms retain counts, proportions, tax, scale
and metadata, including legitimate nested branches and other scientific fields.
Validation checks these objects without changing their values.

Run the commands below from the repository root. Replace `STAGING` with the
staging directory, `ASSETS` with a new empty output directory, and `DATA_TAG` with
a versioned tag such as `data-v1.1.0`. The staged data version must match the tag.
Packaging requires the `zip` command.

```sh
Rscript scripts/validate-study.R STAGING/vandeputte2017
Rscript scripts/package-study.R STAGING/vandeputte2017 .build/one-asset DATA_TAG
Rscript scripts/validate-asset.R .build/one-asset/vandeputte2017.zip STAGING/vandeputte2017
Rscript scripts/build-release.R --input STAGING --tag DATA_TAG --output ASSETS
```

The first three commands allow a single-study check. The final command validates
all staged studies, creates one ZIP per study and writes `studies.candidate.json`
in `ASSETS`. Each ZIP is checked for valid members, checksums, exact extracted
bytes and scientific equality with staging. Existing output files are never
overwritten. The candidate manifest uses status `staged` because local packaging
does not make a study publicly downloadable.

## Uploading and verifying a draft

Before creating a new draft, push the reviewed source commits and authenticate
`gh` with release write access to `Silverman-Lab/mutt`. Then run:

```sh
Rscript scripts/publish-draft.R STAGING ASSETS
```

The script first requires a successful authenticated release listing. It creates
a draft if necessary, uploads missing assets and downloads every asset for
verification. Each download must match the recorded size and checksum, pass
archive checks and reproduce the staged scientific objects, including all named
forms. Only a fully successful run writes `ASSETS/draft-verification.json`.

The same command can resume a draft. It verifies existing assets against staging
and never replaces them. It refuses a published release. A failed authentication
or transfer does not complete verification.

## Publishing and checking public downloads

For `data-v1.0.0`, packaging and draft upload/download verification are already
complete. Retain the existing staging and asset directories for the final checks;
there is no need to rebuild or upload the 40 ZIPs again.

Publication is an explicit maintainer action after review of the draft and its
verification record. For the current draft:

```sh
gh release edit data-v1.0.0 --repo Silverman-Lab/mutt --draft=false
```

After publication, replace `STAGING` and `ASSETS` below with the existing release
directories and verify every pinned URL without authentication:

```sh
Rscript scripts/verify-public.R STAGING ASSETS .build/studies.available.json
```

This script downloads and validates every asset again, then writes a manifest
with the verified studies marked `available`. It also retains existing catalog
entries that are absent from the release, including unavailable studies. Review
the generated file and use it to update `inst/extdata/studies.json`; check that
the unavailable records remain intact before committing the package update.
Never mark a study available before these public checks pass.

Published assets are immutable. Corrections require a new `data-vX.Y.Z` tag and
an updated pinned package manifest. Do not overwrite a published asset or use a
moving `latest` URL.

## Checking the package and release tools

Run the package fixture tests and source build/check when preparing a change for
review. Routine CI checks package code and tiny local fixtures; it must not run
migration, rebuild studies or download study data or FASTQs. Optional
Bioconductor converters are not dependencies of the ordinary download API.

The release-tool fixture test is maintainer-only and lives outside the installed
package:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE); source("tests/testthat/helper-release.R"); testthat::test_file("scripts/tests/test-release-tools.R", stop_on_failure=TRUE)'
```

It packages only tiny local fixtures. Installed-package tests must not depend on
maintainer scripts excluded from the source package.

## Historical v1 capture

The one-time migration established exact agreement with API outputs from commit
`01c648b1f736b1af019d1c0f0231b73de2c21e49`. It did not rebuild taxonomic profiles.
The old package was installed in an isolated R library and read its stored
non-functional study inputs. Original inputs remain in a gitignored recovery
snapshot and unchanged Git history.

These commands document that capture rather than a task still to perform:

```sh
TZ=UTC Rscript scripts/migration/capture-current.R OLD_R_LIBRARY LEGACY_STUDIES PARSER_INDEX .build/capture
Rscript scripts/migration/stage-current.R .build/capture/default .build/capture/aligned .build/data-v1.0.0
```

UTC is required for reproducible capture. The capture script guards DADA2
classification and denoising entry points so that missing processed inputs cause
an explicit failure rather than new scientific preprocessing. Some parsers still
need the `dada2` namespace to read stored tables; that parsing dependency does not
permit classification or denoising. Capture strips only the deliberately removed
empty `function` field and preserves every other scientific field.

Comparison uses exact scientific-object equality, not a numerical tolerance or
a check of dimensions alone. Existing parser failures remain recorded rather
than receiving substituted historical objects. The separate 37-study historical
RData is a recovery reference, not the selected baseline; `stage-v1.R` is retained
only as its historical extraction record. The bounded Dreier scale recovery and
its exact baseline verification are documented in [the migration result](migration/FINAL.md).
The scripts in `scripts/migration/` are historical tools and are not supported as
the normal path for future releases.
