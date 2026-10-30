# Maintaining study data releases

Git stores the R package, study descriptions and processing scripts. GitHub
Releases store one ZIP per processed study, so users download only what they
need. Each ZIP contains saved R objects and a `provenance.json` file describing
where they came from and how their files were checked.

The catalog, `inst/extdata/studies.json`, connects study IDs to data versions,
download URLs, file sizes and checksums. A checksum lets the package check that
a downloaded file is exactly the expected file. Data versions such as
`data-v1.1.0` are independent of R package versions.

## Prepare study files

Prepare new or updated studies with their own processing scripts. Keep temporary
files in the ignored `.build/` directory. Do not rerun scientific processing
just to package existing data or test downloads.

Each prepared folder is named by study ID and contains:

- `study.rds`: the default saved study object;
- `provenance.json`: its study ID, data and format versions, source information,
  preparation decisions and file checksums;
- any additional saved versions declared in the file's `variants` entries,
  with their filenames and checksums.

Preserve scientific values, names, units, missing values and separate results
from different methods. When including another saved version, explain its rules
for matching samples. The name `aligned.rds` alone does not guarantee that every
table contains the same samples in the same order.

## Check and package

Run these commands from the repository root. Replace `STAGING` with the prepared
study folders, `ASSETS` with a new empty output folder and `DATA_TAG` with a
version such as `data-v1.1.0`. Each study's `provenance.json` must record the
corresponding `data_version`, such as `1.1.0`. Packaging requires the `zip` command.

```sh
Rscript scripts/validate-study.R STAGING/vandeputte2017
Rscript scripts/package-study.R STAGING/vandeputte2017 .build/one-study DATA_TAG
Rscript scripts/validate-asset.R .build/one-study/vandeputte2017.zip STAGING/vandeputte2017
Rscript scripts/build-release.R --input STAGING --tag DATA_TAG --output ASSETS
```

The first three commands check one study. The last command checks and packages
all prepared studies and writes `ASSETS/studies.candidate.json`, a proposed
catalog. Checks cover archive contents, file checksums, extracted bytes and
study values. Existing output files are never overwritten. Newly packaged
studies have status `staged` until public downloads have been verified.

## Upload and check a draft

Push the reviewed source commits first. Authenticate `gh` with permission to
write releases in `Silverman-Lab/mutt`, then run:

```sh
Rscript scripts/publish-draft.R STAGING ASSETS
```

This uploads missing ZIPs to a draft release and downloads every ZIP again to
check it against the prepared study. It writes `ASSETS/draft-verification.json`
only when all checks pass. The same command can resume an interrupted upload:
it checks existing ZIPs rather than replacing them, and refuses a published
release.

Keep the small proposed catalog and verification record. Once the uploaded
files have been verified, local bulk copies can be removed. To recover them
later without scientific processing, download the ZIPs into the folder with the
retained `studies.candidate.json`, then restore into a new folder:

```sh
gh release download DATA_TAG --repo Silverman-Lab/mutt --pattern '*.zip' --dir ASSETS
Rscript scripts/restore-staging.R ASSETS NEW_STAGING
```

Restoration checks the downloads and will not overwrite an existing study folder.

## Publish and check public downloads

A maintainer publishes the draft explicitly after review. After publication,
check the public downloads without authentication:

```sh
Rscript scripts/verify-public.R STAGING ASSETS .build/studies.available.json
```

This checks every released study and writes a catalog marking verified studies
`available`. It retains other catalog entries, including unavailable studies.
Review the result, use it to update `inst/extdata/studies.json` and commit the
package update. Never mark studies available before these public checks pass.

Never replace published ZIPs or use a moving `latest` URL. Corrections require
a new data version and an updated package catalog.

## Check the package and release tools

Run the package tests and an R source-package build and check before review.
Routine CI uses small local examples; it must not download study data or run
scientific processing. Optional converters do not add requirements for ordinary
study loading.

Check the maintenance tools separately:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE); source("tests/testthat/helper-release.R"); testthat::test_dir("scripts/tests", stop_on_failure=TRUE)'
```

These tests package and restore small local examples. Installed-package tests
must not depend on maintenance scripts excluded from the source package.
