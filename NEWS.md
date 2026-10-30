# mutt 0.1.0.9001

* Study data are downloaded separately from the R package. Use `mutt_studies()`
  to find studies and `mutt(studies)` to load the ones you need.
* Downloads are checked before use and saved on your computer for reuse. Use
  `refresh = TRUE` to download again, or `mutt_cache_remove()` to remove a copy.
* Study objects retain their counts, relative abundances, taxonomy, load
  measurements, sample information and separate results from different methods.
  `representation = "aligned"` loads a second saved version using each study's
  existing rules for matching samples.
* `mutt_source()`, `mutt_recipe()` and `mutt_download_raw()` provide source
  information, processing-script links and original inputs separately.
* Optional phyloseq and miaverse converters accept a selected counts table and
  match its sample and feature IDs to the supplied information tables.
* The empty `function` field and functional-inference implementation are removed.
