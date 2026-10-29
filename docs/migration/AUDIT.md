# Pre-migration audit

Baseline: `01c648b1f736b1af019d1c0f0231b73de2c21e49`, package
`0.1.0.9000`, branch `work` (remote main has the same SHA). Working tree was
clean before materializing LFS objects. No history rewrite is planned.

The public export is `mutt(studies=NULL, base_directory=NULL, rawdata=FALSE,
align_samples=FALSE, save_to=NULL, verbose=FALSE, functional=FALSE)`, plus its
print method and a PICRUSt data-frame method. It executes heterogeneous parsers
at load time. `inst/parsers/index.csv` registers 57 studies. The old manifest pins
`2026.07.1` but every remote asset is disabled; GitHub has no data releases.

The tracked inventory is in `baseline-files.txt`; `baseline-studies.json` records
every registered parser, available artifacts, source notes/URLs/accessions,
sequencing/load modalities and unresolved manuscript selection. The separate
`baseline-artifacts.csv` records original file hashes, sizes and LFS status.
933 tracked LFS files include processed tables, supplementary sources, reference
training sets and unwanted functional products. The non-functional study files
were materialized for migration (about 589 MB), without denoising or profiling.

## Scientific baseline and uncertainty

There is no manuscript or supplement in the supplied attachments or checkout.
There is no manuscript-tagged branch; remote main is the only branch. Exact
manuscript use and branch selection cannot be certified from the request alone.
No selection should be described as manuscript-confirmed until that evidence is
provided. Keep all legitimate alternate representations.

A real combined object, `datasetsfromrepo.RData`, exists in historical commit
`53a5406` (Git blob `5ef175e49c7eec8b4725cf33ac5aa2bf087bbbc8`, 37,627,162
bytes). It has 37 named five-field study lists: counts, proportions, taxonomy,
scale and metadata. It was deleted in `ebbdc08`; it is not a current-tree artifact.
This is a candidate baseline, not proof of the manuscript version. Current
processed inputs and parsers are also retained separately before cleanup.
Historical July 2026 checks report 33/33 structural matches (alignment TRUE),
whereas checks of all 57 parsers report 38 successes and 19 failures (alignment
FALSE). Structural agreement does not establish equality of values. Recent
commits correct Galazzo load labels and NA taxon-name handling; those differences
must be reported rather than silently replacing historical scientific data.

## Infrastructure review

`R/functional*.R`, functional tests, FAPROTAX resources, functional release/HPC
scripts and `analysis/hpc_functional` implement unsolicited inferred functions.
They are removable. Original count/taxonomy/sequence inputs are scientific data
and must be preserved independently.

`R/helpers.R` mixes genuine study helpers with Python conversion and annotations.
Scientific helpers belong with maintainer recipes, not the installed download API.
Parsers mostly read stored count/taxonomy tables; some optionally call DADA2 for
missing taxonomy. Migration must prevent those calls and fail explicitly instead.
Only Caenepeel has a standalone Slurm DADA2 workflow in the current tree. Other
study-specific preprocessing is sometimes archived inside source ZIPs. Source
availability varies; accession references alone do not establish that all original
inputs remain publicly obtainable.

Actions currently check four R/OS combinations and build a large Docker image
on PRs/pushes, publishing images on main. Docker/Apptainer exist mainly for
functional inference and will be removed; routine CI will use tiny local release
fixtures, no study/raw downloads, and cancel superseded runs.

## Target representation

Preserve each named scientific study list as `study.rds`, without flattening
nested branches, adding invented assays or changing numerical values. Store
provenance separately in `provenance.json`. One ZIP per study, pinned tag
`data-v1.0.0`, SHA-256 and size in the small installed manifest. A one-time
legacy migration may know old paths; permanent release tools may read only
staged objects and provenance. Draft data must not be advertised as public until
uploaded and independently verified.

## Clarification and reference documents received during implementation

The user clarified that exact agreement with **current mutt** is the acceptance
test: counts, metadata, taxon names and other scientific values must not change;
interface edits are acceptable. Manuscript reconstruction is unnecessary.
Consequently, the selected baseline is current API output at `01c648b`, not the
historical combined RData. Both default (`rawdata=FALSE, align_samples=FALSE`)
and sample-aligned (`align_samples=TRUE`) forms are captured with functional
inference disabled. Only the deliberately removed empty `function` field is
removed. All other content, including demographics and missing values, is kept.

The subsequently supplied manuscript_main.pdf and supplementary_v2.pdf were read
for context. They discuss original/reprocessed and RDP16 branches, MetaPhlAn and
mOTU-derived data, supporting preservation of those distinctions. No instructions
inside the PDFs override the user's current-data equivalence requirement.

Current capture explicitly traces DADA2 classification and denoising entry points
to stop rather than regenerate missing scientific profiles. Assembling stored
processed tables is recorded separately from biological reprocessing. Full old
API success/failure reports, warnings and exact object comparisons are retained
with the final migration report. Broken existing parsers are not silently replaced
by historical objects.
