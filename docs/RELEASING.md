# Releasing study data

Use this guide when adding or correcting study data in MUTT. You will package
prepared R data, upload them to a GitHub release, then update the package's
study catalog so users can download them. Changes to package code or help alone
do not require a new data release.

The examples use `data-v1.1.0`. Choose a new version for your release: published
study files stay available under their original versions, so corrections go in
a new release rather than replacing files used by earlier analyses.

Run the commands from the repository root. You need R with `digest` and
`jsonlite`, the `zip` command, and GitHub CLI (`gh`) authenticated with permission
to create releases in `Silverman-Lab/mutt`.

## 1. Collect the prepared studies

Start with the completed study files produced by the processing scripts in
`studies/`. The release tools bundle those saved files for download.

For this example, put the studies you want to release under
`.build/data-v1.1.0/prepared/`, with one folder per study:

```text
.build/data-v1.1.0/prepared/
  vandeputte2017/
    study.rds
    provenance.json
  contijoch2019/
    study.rds
    provenance.json
```

`study.rds` holds the prepared study tables. `provenance.json` records their
sources, preparation and file checksums. The processing scripts create both
files. Before packaging, set `data_version` in `provenance.json` to `1.1.0` to
match the example release tag `data-v1.1.0`. If a study has additional saved RDS
files listed in `provenance.json`, include those too.

Only include studies being released in this version; the catalog update below
retains entries for other studies. The `.build/` folder is excluded from Git,
so these working copies and the ZIPs created next will not be committed with
package code. If you need copies from an existing release, see recovery below.

## 2. Build the study ZIPs

```sh
Rscript scripts/build-release.R \
  --input .build/data-v1.1.0/prepared \
  --tag data-v1.1.0 \
  --output .build/data-v1.1.0/release
```

This checks each study, creates one ZIP per study and checks that extracting it
reproduces the input files and study values. The output folder also contains
`studies.candidate.json`: the proposed catalog entries with download addresses,
file sizes and checksums. Its entries have status `staged`, meaning the data
have been prepared but are not yet ready for public downloads.

The output folder must be new or empty. If the build fails, correct the reported
problem and retry with a new or empty output folder. You can check a particular
input on its own with `Rscript scripts/validate-study.R PATH_TO_STUDY_FOLDER`.

## 3. Upload a draft for review

Commit and push the source changes that belong with the release before uploading.
The upload script uses the current commit when creating the draft:

```sh
Rscript scripts/publish-draft.R \
  .build/data-v1.1.0/prepared .build/data-v1.1.0/release
```

The script uploads missing ZIPs, then downloads each one and compares it with
your prepared study. Once every check passes, it writes
`.build/data-v1.1.0/release/draft-verification.json`. The GitHub release remains
a draft for review.

If a transfer fails, rerun the same command. It checks files already uploaded
and adds missing files. A comparison failure needs investigation before you
continue: the uploaded file differs from the prepared copy. The script will
not replace uploaded files or modify a published release.

## 4. Publish and enable downloads in the package

Review the draft's study files and release notes in GitHub Releases, then publish
it there. Next, check that users can download the files without signing in:

```sh
Rscript scripts/verify-public.R \
  .build/data-v1.1.0/prepared .build/data-v1.1.0/release \
  .build/data-v1.1.0/studies.available.json
```

A successful run writes `studies.available.json`, with the checked studies
marked `available` and the other catalog entries retained. Review that file,
replace `inst/extdata/studies.json` with it, and commit the package update. Users
receive the new download information when they update the package.

If a public download fails, the script stops. Keep the existing package catalog
until the problem is resolved and the command succeeds.

## Recover local copies from a release

The upload and public-download checks compare the released files with your
prepared studies, so keep both local folders until those checks are complete.
Afterward, retain `studies.candidate.json` and `draft-verification.json`; the
large ZIPs and prepared copies can be removed to free disk space.

If you later need those copies again, download the ZIPs into the folder containing
that release's `studies.candidate.json`, then restore the study folders:

```sh
gh release download data-v1.1.0 --repo Silverman-Lab/mutt \
  --pattern '*.zip' --dir .build/data-v1.1.0/release
Rscript scripts/restore-staging.R \
  .build/data-v1.1.0/release .build/data-v1.1.0/restored
```

Restoration checks each download against the saved catalog. It writes the study
folders to `restored/`, which must not already exist. Use that folder in place
of `prepared/` for any further release checks.

## Check changes to the release tools

When editing the release scripts, run their tests with `pkgload`, `testthat`,
`withr`, `digest` and `jsonlite` installed, and the `zip` command available:

```sh
Rscript -e 'pkgload::load_all(".", quiet=TRUE); source("tests/testthat/helper-release.R"); testthat::test_dir("scripts/tests", stop_on_failure=TRUE)'
```

These tests build and restore small example studies on your computer. They do
not upload files or publish a release.
