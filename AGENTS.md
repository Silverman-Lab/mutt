# Working on *mutt*

*mutt* provides processed microbiome study data from R. Users list studies with
`mutt_studies()` and load the ones they need with `mutt(studies)`.

Write the package name as italic lowercase *mutt* in prose. Keep function names,
paths and other code identifiers unchanged.

## Project layout

- `R/`: installed package code for finding, downloading and loading studies.
- `studies/`: source descriptions and study-specific download and processing scripts.
- `scripts/`: reusable maintenance and release tools.
- `tests/`: small examples for checking package behavior.
- [docs/RELEASING.md](docs/RELEASING.md): instructions for maintaining data releases.

GitHub Releases hold the study ZIPs. `inst/extdata/studies.json` records their
versions, URLs, sizes and checksums. Keep large data outside the active source
tree. Use the ignored `.build/` directory for temporary work and clean it up
when finished.

## Preserve study data

Keep scientific values, names, units, missing values and separate results from
different methods intact. Downloading and packaging must not change counts or
load measurements, or choose a table implicitly.

The default saved study may have different samples in different tables. The
optional `aligned` version follows each study's existing sample-ID matching
rules. Loading reads a saved object; it does not rerun processing or guarantee
that all tables contain the same samples in the same order.

Keep scientific processing study-specific. Shared helpers handle retrieval and
file checks. Describe missing inputs honestly; never invent sources or results.
Preserve original source notes and descriptions of the authors' methods.

## Downloads and releases

Use explicit data versions rather than a moving `latest` URL. Check download
sizes, checksums, archive contents and stored files before use. Extract into a
temporary directory before moving a verified study into the user's saved-data
folder. A failed refresh must leave the previous copy usable.

Release tools package prepared objects without running scientific processing.
Never overwrite published files: corrections require a new data version. Mark
studies available only after publication and successful public-download checks.

## Dependencies and checks

Keep installed dependencies small. Processing dependencies belong in study
scripts. The installed package must work independently of maintenance scripts.
Tests use small local examples; routine CI must not download or rebuild studies.
Functional inference and downstream biological analysis are outside *mutt*'s scope.
Never rewrite Git history as part of cleanup.

Write API help in the roxygen comments in `R/` and regenerate `man/` with
`roxygen2::roxygenise(roclets = "rd")`. Do not edit generated help separately.
