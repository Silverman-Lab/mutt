# mutt 0.1.0.9000

* Added `functional = "REVALIDATE"` to reconstruct validated PICRUSt2 results
  from retained normal or `.failed/raw` output without launching PICRUSt2.
  Pathway-only sample omissions are zero-filled and reported as
  `generated_with_warning`; EC/KO disagreement and unexpected pathway samples
  remain failures.
* The HPC study runner now accepts `revalidate` and exits with status 11 when
  any functional manifest branch has `status = "failed"`, preventing failed
  studies from receiving completion checkpoints.
* Functional inference now routes only explicitly identified amplicon branches
  to PICRUSt2 and FAPROTAX; shotgun, metagenomic, and unknown-modality branches
  are retained in the manifest as skipped.
* PICRUSt2 EC and KO outputs must agree exactly. Pathway outputs may omit EC/KO
  samples, which are restored as zero-prediction rows with explicit
  reconciliation QC; unexpected pathway samples remain validation failures.
* DADA2 PICRUSt2 inputs now merge exact reverse-complement ASV duplicates before
  inference while preserving every sample's read total and rebuilding aligned
  consensus taxonomy for the merged features.
* Added safe, explicit one-file study-data archiving with Git LFS verification.
* Study TSV and RDS readers now support either an uncompressed file or its
  one-file `.zip` archive, preserving local rebuild workflows and archived Git
  checkouts.
* Functional caches can be restored from validated `functional.zip` archives;
  cached PICRUSt2 output paths are rebased to the current checkout.
* Successful PICRUSt2 caches discard unneeded orientation and intermediate
  alignment directories after parsing while retaining predictions, stratified
  contribution files, mappings, taxonomy, QC, logs, and provenance.
* Main-branch container builds are published to GitHub Container Registry for
  conversion and execution with Apptainer on HPC systems.
* Restructured MUTT as an installable package with `mutt()` as its sole public
  API, bundled parsers, local/remote study resolution, structured audit records,
  and checksum-verified atomic caching.
* Added all-study data-release readiness, build, local verification, and remote
  verification tooling. Production remote flags remain disabled until assets
  are published and verified.
* Added optional cached PICRUSt2 and FAPROTAX inference through
  `mutt(functional = FALSE)`, with `TRUE` and `"REBUILD"` modes.
* Functional files are stored under each study's `functional/` directory and
  returned in the study's `function` entry with branch-level provenance.
* PICRUSt2 inputs now undergo per-ASV strand selection against its bacterial
  and archaeal reference HMMs, with alignment fractions and decisions retained.
* PICRUSt2 now always requests stratified EC, KO, and MetaCyc output. Large
  ASV-contribution tables remain file-backed and are accessed through
  `as.data.frame()` on a returned PICRUSt2 branch.
* Bundled FAPROTAX 1.2.12 files are discovered without directory arguments.
* FAPROTAX now prefers matched count and taxonomy branches and falls back to
  matched proportion and taxonomy branches when counts are unavailable.
* Preserved valid sequence columns when taxonomy row names are ASV identifiers.
* Added a compatibility audit for the 33 study outputs present in the saved
  pre-restructure validation artifact. The audit records that this baseline
  used sample alignment and keeps it separate from the default-mode parser
  inventory.

# mutt 0.1.0

* Initial CRAN submission.
