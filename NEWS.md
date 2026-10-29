# mutt (development)

* Replace parser-on-load use with a small pinned-release download API, user cache,
  checksum validation, explicit refresh/removal and source/recipe discovery.
* Preserve existing scientific objects and sample-aligned forms without rerunning
  taxonomic profiling. Failed legacy parsers remain explicitly unavailable.
* Separate one-time legacy migration from permanent staged release tooling.
* Standardize study recipe interfaces while preserving scientific internals.
* Remove downstream functional inference and its containers/HPC infrastructure.
* Add optional phyloseq and miaverse conversion of explicitly selected tables.
