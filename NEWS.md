# mutt 0.1.0.9001

* MUTT now separates the R package from processed study data. `mutt_studies()`
  lists the catalog, and `mutt(studies)` retrieves only requested studies from a
  pinned data release. Downloads are checked with SHA-256 and stored in a local
  cache, with explicit refresh and removal options.
* Study objects preserve the scientific outputs of MUTT at commit `01c648b`,
  including nested branches, metadata and missing values, without rerunning
  taxonomic profiling. The default form matches the old API's `rawdata=FALSE,
  align_samples=FALSE` output; `representation="aligned"` selects its stored
  `align_samples=TRUE` form.
* The catalog records 40 usable processed studies, 17 existing parser failures
  and six additional source records without processed outputs. All 40 study ZIPs
  have passed upload/download verification in the `data-v1.0.0` draft. The draft
  remains unpublished; the package will mark these studies available only after
  publication and public URL verification.
* Source discovery and study recipes now have a common interface while retaining
  study-specific processing. Maintainer release tools validate and package
  staged objects separately from the one-time historical migration.
* The empty `function` field and downstream functional-inference implementation,
  resources, containers and HPC infrastructure have been removed. Inferred
  functional profiles and pathways remain outside MUTT's scope.
* Optional phyloseq and miaverse converters work with explicitly selected tables
  and check sample and taxon alignment before conversion.
