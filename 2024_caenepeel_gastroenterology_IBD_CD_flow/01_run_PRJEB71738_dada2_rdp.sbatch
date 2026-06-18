#!/usr/bin/env bash
#SBATCH --job-name=PRJEB71738_dada2_rdp
#SBATCH --account=one_sc_default
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=PRJEB71738_dada2_rdp_%j.out
#SBATCH --error=PRJEB71738_dada2_rdp_%j.err

set -euo pipefail

############################################
# User-editable settings
############################################

ACC="${ACC:-PRJEB71738}"
ENV_NAME="${ENV_NAME:-dada2_rdp}"

PROJECT_DIR="${PROJECT_DIR:-$PWD/${ACC}_dada2_rdp}"
RAW_DIR="${PROJECT_DIR}/raw_fastq"
TRIM_DIR="${PROJECT_DIR}/trimmed_fastq"
FILTER_DIR="${PROJECT_DIR}/filtered_fastq"
META_DIR="${PROJECT_DIR}/metadata"
REF_DIR="${PROJECT_DIR}/refs"
SCRIPT_DIR="${PROJECT_DIR}/scripts"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"

THREADS="${SLURM_CPUS_PER_TASK:-24}"

# Starting values for V4 16S MiSeq 2x250 after primer removal.
# Tune these if read_tracking.tsv shows poor filtering or poor merging.
# After cutadapt removes the 515F/806R primers from 2x250 V4 reads,
# the usable reads are typically ~230 bp, not 250 bp. The previous
# TRUNC_F=240 could therefore discard every forward read.
TRUNC_F="${TRUNC_F:-220}"
TRUNC_R="${TRUNC_R:-180}"
MAXEE_F="${MAXEE_F:-2}"
MAXEE_R="${MAXEE_R:-5}"
MIN_OVERLAP="${MIN_OVERLAP:-12}"

# 515F / 806R primers. Script handles either orientation.
FWD_PRIMER="${FWD_PRIMER:-GTGYCAGCMGCCGCGGTAA}"
REV_PRIMER="${REV_PRIMER:-GGACTACNVGGGTWTCTAAT}"
FWD_RC="${FWD_RC:-TTACCGCGGCKGCTGRCAC}"
REV_RC="${REV_RC:-ATTAGAWACCCBNGTAGTCC}"

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$FILTER_DIR" "$META_DIR" "$REF_DIR" "$SCRIPT_DIR" "$RESULTS_DIR" "$LOG_DIR"

echo "[job] Started: $(date)"
echo "[job] Host: $(hostname)"
echo "[job] Accession: $ACC"
echo "[job] Project dir: $PROJECT_DIR"
echo "[job] Threads: $THREADS"
echo "[job] Memory requested: ${SLURM_MEM_PER_NODE:-unknown}"
echo "[job] CPUs requested: ${SLURM_CPUS_PER_TASK:-unknown}"

############################################
# Activate conda
############################################

if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
else
  echo "ERROR: conda not found. Load your HPC conda module before sbatch, or edit this script." >&2
  exit 1
fi

conda activate "$ENV_NAME"

echo "[job] Conda env:"
conda info --envs | sed -n '1,20p'

echo "[job] Tool versions:"
Rscript -e 'cat("R:", R.version.string, "\n"); cat("dada2:", as.character(packageVersion("dada2")), "\n")'
cutadapt --version || true

############################################
# Make sure references exist
############################################

download_if_missing() {
  local url="$1"
  local out="$2"

  if [[ -s "$out" ]]; then
    echo "[refs] Already exists: $out"
  else
    echo "[refs] Downloading: $out"
    curl -L --fail --retry 5 --retry-delay 10 -o "$out" "$url"
  fi
}

download_if_missing \
  "https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz?download=1" \
  "${REF_DIR}/rdp_train_set_16.fa.gz"

download_if_missing \
  "https://zenodo.org/records/14168771/files/rdp_19_toGenus_trainset.fa.gz?download=1" \
  "${REF_DIR}/rdp_19_toGenus_trainset.fa.gz"

############################################
# Download ENA metadata
############################################

echo "[metadata] Downloading ENA read_run file report"

META="${META_DIR}/${ACC}_ena_read_run.tsv"
FIELDS="run_accession,sample_accession,secondary_sample_accession,experiment_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout,library_strategy,library_source,instrument_model,read_count,base_count"

curl -L --fail --retry 5 \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=${FIELDS}&format=tsv" \
  -o "$META"

if [[ "$(wc -l < "$META")" -le 1 ]]; then
  echo "ERROR: ENA returned no FASTQ metadata for accession: $ACC" >&2
  exit 1
fi

echo "[metadata] Metadata rows including header: $(wc -l < "$META")"

############################################
# Build manifest and download FASTQs
############################################

echo "[download] Building FASTQ URL list and manifest"

python - "$META" "$RAW_DIR" "${META_DIR}/fastq_urls.txt" "${META_DIR}/fastq_md5.txt" "${META_DIR}/raw_manifest.tsv" <<'PY'
import csv
import os
import sys
from itertools import zip_longest

meta_file, raw_dir, urls_file, md5_file, manifest_file = sys.argv[1:6]

seen = set()

with open(meta_file, newline="") as f, \
     open(urls_file, "w") as urls_out, \
     open(md5_file, "w") as md5_out, \
     open(manifest_file, "w", newline="") as man_out:

    reader = csv.DictReader(f, delimiter="\t")

    writer = csv.writer(man_out, delimiter="\t")
    writer.writerow([
        "sample_id",
        "run_accession",
        "sample_accession",
        "secondary_sample_accession",
        "experiment_accession",
        "layout",
        "forward",
        "reverse"
    ])

    for row in reader:
        run = row.get("run_accession", "").strip()
        sample = row.get("sample_accession", "").strip()
        secondary_sample = row.get("secondary_sample_accession", "").strip()
        experiment = row.get("experiment_accession", "").strip()
        layout = row.get("library_layout", "").strip().upper()

        urls = [x.strip() for x in row.get("fastq_ftp", "").split(";") if x.strip()]
        md5s = [x.strip() for x in row.get("fastq_md5", "").split(";") if x.strip()]

        if not urls:
            continue

        if layout != "PAIRED" or len(urls) != 2:
            raise SystemExit(
                f"Run {run} is not paired-end with exactly 2 FASTQ files. "
                f"layout={layout}, n_fastq={len(urls)}"
            )

        local_files = []
        for url, md5 in zip_longest(urls, md5s, fillvalue=""):
            full_url = url if url.startswith(("http://", "https://", "ftp://")) else "https://" + url
            local = os.path.abspath(os.path.join(raw_dir, os.path.basename(url)))

            if full_url not in seen:
                urls_out.write(full_url + "\n")
                seen.add(full_url)

            local_files.append(local)

            if md5:
                md5_out.write(f"{md5}  {local}\n")

        writer.writerow([
            run,
            run,
            sample,
            secondary_sample,
            experiment,
            "PAIRED",
            local_files[0],
            local_files[1]
        ])
PY

echo "[download] Downloading FASTQs with parallel curl"
echo "[download] Number of FASTQ files: $(wc -l < "${META_DIR}/fastq_urls.txt")"

cat "${META_DIR}/fastq_urls.txt" | \
  xargs -n 1 -P 6 -I {} bash -c '
    url="$1"
    out_dir="$2"
    out="${out_dir}/$(basename "$url")"
    if [[ -s "$out" ]]; then
      echo "[download] Already exists: $out"
    else
      echo "[download] Downloading: $url"
      curl -L --fail --retry 5 --retry-delay 10 --continue-at - -o "$out" "$url"
    fi
  ' _ {} "$RAW_DIR"

if [[ -s "${META_DIR}/fastq_md5.txt" ]]; then
  echo "[download] Checking MD5 sums, but not stopping on checksum failures"
  set +e
  md5sum -c "${META_DIR}/fastq_md5.txt" > "${LOG_DIR}/md5sum_check.log" 2>&1
  MD5_STATUS=$?
  set -e

  cat "${LOG_DIR}/md5sum_check.log"

  if [[ "$MD5_STATUS" -ne 0 ]]; then
    echo "[download] WARNING: one or more MD5 checksums failed."
    echo "[download] Continuing anyway because checksum failure is non-fatal in this version."
    grep "FAILED" "${LOG_DIR}/md5sum_check.log" || true
  fi
fi

############################################
# Validate FASTQs and drop corrupt pairs
############################################

echo "[validate] Checking FASTQ integrity with gzip -t before cutadapt"

SKIPPED_LOG="${LOG_DIR}/skipped_corrupt_fastq.log"
: > "$SKIPPED_LOG"   # truncate / create

python - "${META_DIR}/raw_manifest.tsv" "${META_DIR}/raw_manifest_clean.tsv" "$SKIPPED_LOG" <<'PY'
import csv
import subprocess
import sys

manifest_in, manifest_out, skipped_log = sys.argv[1:4]

def gzip_ok(path):
    """Return True if the file exists, is non-empty, and passes gzip integrity check."""
    import os
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    result = subprocess.run(
        ["gzip", "-t", path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return result.returncode == 0

kept = 0
skipped = 0

with open(manifest_in, newline="") as fin, \
     open(manifest_out, "w", newline="") as fout, \
     open(skipped_log, "a") as skip_f:

    reader = csv.DictReader(fin, delimiter="\t")
    writer = csv.DictWriter(fout, delimiter="\t", fieldnames=reader.fieldnames)
    writer.writeheader()

    for row in reader:
        fwd = row["forward"]
        rev = row["reverse"]
        run = row["run_accession"]

        fwd_ok = gzip_ok(fwd)
        rev_ok = gzip_ok(rev)

        if fwd_ok and rev_ok:
            writer.writerow(row)
            kept += 1
        else:
            bad = []
            if not fwd_ok:
                bad.append(fwd)
            if not rev_ok:
                bad.append(rev)
            msg = f"[validate] SKIPPING run {run}: corrupt or missing file(s): {', '.join(bad)}"
            print(msg, flush=True)
            skip_f.write(msg + "\n")
            skipped += 1

print(f"[validate] Kept {kept} run(s), skipped {skipped} run(s).", flush=True)

if kept == 0:
    raise SystemExit("ERROR: No valid FASTQ pairs remain after integrity check. Aborting.")
PY

echo "[validate] Skipped runs logged to: $SKIPPED_LOG"

# Use the clean manifest for all downstream steps
RAW_MANIFEST_FOR_TRIM="${META_DIR}/raw_manifest_clean.tsv"

############################################
# Primer trimming with cutadapt
############################################

echo "[cutadapt] Writing cutadapt command file"

python - "${RAW_MANIFEST_FOR_TRIM}" "$TRIM_DIR" "${META_DIR}/trimmed_manifest.tsv" "${SCRIPT_DIR}/cutadapt_commands.sh" \
  "$FWD_PRIMER" "$REV_PRIMER" "$FWD_RC" "$REV_RC" "$THREADS" <<'PY'
import csv
import os
import shlex
import sys

manifest, trim_dir, trimmed_manifest, cmd_file = sys.argv[1:5]
fwd, rev, fwd_rc, rev_rc, threads = sys.argv[5:10]

with open(manifest, newline="") as f, \
     open(trimmed_manifest, "w", newline="") as out, \
     open(cmd_file, "w") as cmds:

    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames

    writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    cmds.write("#!/usr/bin/env bash\n")
    cmds.write("set -euo pipefail\n\n")

    for row in reader:
        sid = row["sample_id"]

        out_f = os.path.abspath(os.path.join(trim_dir, f"{sid}_R1.trim.fastq.gz"))
        out_r = os.path.abspath(os.path.join(trim_dir, f"{sid}_R2.trim.fastq.gz"))

        # Handles normal and reversed R1/R2 orientation.
        # Removes 5-prime primers and possible 3-prime read-through primer complements.
        cmd = [
            "cutadapt",
            "-j", threads,
            "--error-rate", "0.10",
            "--overlap", "10",
            "--minimum-length", "50",
            "-g", "^" + fwd,
            "-g", "^" + rev,
            "-a", fwd_rc,
            "-a", rev_rc,
            "-G", "^" + rev,
            "-G", "^" + fwd,
            "-A", rev_rc,
            "-A", fwd_rc,
            "-o", out_f,
            "-p", out_r,
            row["forward"],
            row["reverse"]
        ]

        cmds.write(" ".join(shlex.quote(x) for x in cmd) + "\n")

        row["forward"] = out_f
        row["reverse"] = out_r
        writer.writerow(row)
PY

chmod +x "${SCRIPT_DIR}/cutadapt_commands.sh"

echo "[cutadapt] Running cutadapt"
bash "${SCRIPT_DIR}/cutadapt_commands.sh" > "${LOG_DIR}/cutadapt.log" 2>&1

############################################
# Write DADA2 R pipeline
############################################

cat > "${SCRIPT_DIR}/run_dada2_rdp_genus.R" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(dada2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 11) {
  stop("Expected args: manifest results_dir filtered_dir rdp16 rdp19 threads truncF truncR maxEEF maxEER minOverlap")
}

manifest_file <- args[[1]]
results_dir <- args[[2]]
filtered_dir <- args[[3]]
rdp16 <- args[[4]]
rdp19 <- args[[5]]
threads <- as.integer(args[[6]])
truncF <- as.integer(args[[7]])
truncR <- as.integer(args[[8]])
maxEEF <- as.numeric(args[[9]])
maxEER <- as.numeric(args[[10]])
minOverlap <- as.integer(args[[11]])

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)

mt <- if (.Platform$OS.type == "windows") FALSE else threads

set.seed(100)

message("[R] Reading manifest: ", manifest_file)

manifest <- read.delim(
  manifest_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(manifest) == 0) {
  stop("Manifest has zero rows.")
}

sample_names <- make.unique(manifest$sample_id)

fnFs <- manifest$forward
fnRs <- manifest$reverse

names(fnFs) <- sample_names
names(fnRs) <- sample_names

if (!all(file.exists(fnFs))) {
  stop("Some forward FASTQs are missing.")
}

if (!all(file.exists(fnRs))) {
  stop("Some reverse FASTQs are missing.")
}

filtFs <- file.path(filtered_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample_names, "_R_filt.fastq.gz"))

names(filtFs) <- sample_names
names(filtRs) <- sample_names

fastq_lengths_preview <- function(files, n_records = 200, max_files = 20) {
  files <- files[seq_len(min(length(files), max_files))]
  out <- integer(0)
  for (ff in files) {
    con <- gzfile(ff, open = "rt")
    on.exit(close(con), add = TRUE)
    x <- readLines(con, n = n_records * 4, warn = FALSE)
    close(con)
    on.exit(NULL, add = FALSE)
    if (length(x) >= 4) {
      seq_lines <- x[seq(2, length(x), by = 4)]
      out <- c(out, nchar(seq_lines))
    }
  }
  out
}

lenF <- fastq_lengths_preview(fnFs)
lenR <- fastq_lengths_preview(fnRs)
if (length(lenF) > 0 && length(lenR) > 0) {
  message("[R] trimmed read length preview F: min/median/max = ",
          min(lenF), "/", median(lenF), "/", max(lenF))
  message("[R] trimmed read length preview R: min/median/max = ",
          min(lenR), "/", median(lenR), "/", max(lenR))
  if (truncF > 0 && truncF > min(lenF)) {
    message("[R] Requested truncF=", truncF,
            " is longer than at least one previewed forward read. Lowering to ", min(lenF), ".")
    truncF <- min(lenF)
  }
  if (truncR > 0 && truncR > min(lenR)) {
    message("[R] Requested truncR=", truncR,
            " is longer than at least one previewed reverse read. Lowering to ", min(lenR), ".")
    truncR <- min(lenR)
  }
}
message("[R] Using truncLen = c(", truncF, ", ", truncR, ") and maxEE = c(", maxEEF, ", ", maxEER, ")")

message("[R] filterAndTrim")
filt_out <- filterAndTrim(
  fwd = fnFs,
  filt = filtFs,
  rev = fnRs,
  filt.rev = filtRs,
  truncLen = c(truncF, truncR),
  maxN = 0,
  maxEE = c(maxEEF, maxEER),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = mt
)

keep <- filt_out[, "reads.out"] > 0 & file.exists(filtFs) & file.exists(filtRs)
if (!all(keep)) {
  dropped <- data.frame(
    SampleID = sample_names[!keep],
    run_accession = manifest$run_accession[!keep],
    reads_in = filt_out[!keep, "reads.in"],
    reads_out = filt_out[!keep, "reads.out"],
    stringsAsFactors = FALSE
  )
  write.table(
    dropped,
    file = file.path(results_dir, "dropped_zero_after_filtering.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  warning(sum(!keep), " sample(s) had zero reads after filtering and were dropped before denoising.")
}

if (!any(keep)) {
  stop("No samples have reads after filterAndTrim. Check truncLen/maxEE and primer/orientation handling.")
}

manifest <- manifest[keep, , drop = FALSE]
sample_names <- sample_names[keep]
fnFs <- fnFs[keep]
fnRs <- fnRs[keep]
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]
filt_out <- filt_out[keep, , drop = FALSE]

message("[R] Samples retained after filtering: ", length(sample_names))

message("[R] learnErrors forward")
errF <- learnErrors(filtFs, multithread = mt)

message("[R] learnErrors reverse")
errR <- learnErrors(filtRs, multithread = mt)

message("[R] dada forward")
dadaFs <- dada(filtFs, err = errF, multithread = mt)

message("[R] dada reverse")
dadaRs <- dada(filtRs, err = errR, multithread = mt)

getN <- function(x) sum(getUniques(x))

message("[R] mergePairs")
mergers <- mergePairs(
  dadaFs,
  filtFs,
  dadaRs,
  filtRs,
  minOverlap = minOverlap,
  maxMismatch = 0,
  verbose = TRUE
)

merged_n <- sapply(mergers, getN)
if (all(merged_n == 0)) {
  stop("No read pairs merged. This usually means R1/R2 orientation is swapped/reverse-complemented, truncLen leaves too little overlap, or minOverlap is too strict.")
}

message("[R] makeSequenceTable")
seqtab <- makeSequenceTable(mergers)

message("[R] removeBimeraDenovo")
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = mt,
  verbose = TRUE
)

if (ncol(seqtab_nochim) == 0) {
  stop("No ASVs retained after chimera removal.")
}

track <- data.frame(
  SampleID = sample_names,
  run_accession = manifest$run_accession,
  sample_accession = manifest$sample_accession,
  input = filt_out[, "reads.in"],
  filtered = filt_out[, "reads.out"],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab_nochim),
  retained_fraction = rowSums(seqtab_nochim) / filt_out[, "reads.in"],
  stringsAsFactors = FALSE
)

seqs <- colnames(seqtab_nochim)
asv_ids <- paste0("ASV", seq_along(seqs))

seqtab_asv <- seqtab_nochim
colnames(seqtab_asv) <- asv_ids

asv_sequences <- data.frame(
  ASV = asv_ids,
  Sequence = seqs,
  Length = nchar(seqs),
  stringsAsFactors = FALSE
)

counts_by_run <- data.frame(
  SampleID = rownames(seqtab_asv),
  seqtab_asv,
  check.names = FALSE
)

write_tsv <- function(x, file) {
  write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

write_tsv(manifest, file.path(results_dir, "sample_manifest_used.tsv"))
write_tsv(track, file.path(results_dir, "read_tracking.tsv"))
write_tsv(asv_sequences, file.path(results_dir, "asv_sequences.tsv"))
write_tsv(counts_by_run, file.path(results_dir, "asv_counts_by_run.tsv"))

if ("sample_accession" %in% names(manifest) && all(nzchar(manifest$sample_accession))) {
  sample_accessions <- manifest$sample_accession[match(rownames(seqtab_asv), sample_names)]

  counts_by_sample_mat <- rowsum(
    seqtab_asv,
    group = sample_accessions,
    reorder = FALSE
  )

  counts_by_sample <- data.frame(
    SampleID = rownames(counts_by_sample_mat),
    counts_by_sample_mat,
    check.names = FALSE
  )

  write_tsv(counts_by_sample, file.path(results_dir, "asv_counts_by_sample_accession.tsv"))
}

saveRDS(seqtab_nochim, file.path(results_dir, "seqtab_nochim_sequence_columns.rds"))
saveRDS(seqtab_asv, file.path(results_dir, "seqtab_nochim_asv_columns.rds"))
saveRDS(errF, file.path(results_dir, "errF.rds"))
saveRDS(errR, file.path(results_dir, "errR.rds"))

format_tax <- function(tax, label) {
  tax_df <- as.data.frame(unname(tax), stringsAsFactors = FALSE)
  names(tax_df) <- colnames(tax)

  out <- data.frame(
    ASV = asv_ids,
    Sequence = seqs,
    tax_df,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  write_tsv(out, file.path(results_dir, paste0("taxonomy_", label, ".tsv")))
  saveRDS(tax, file.path(results_dir, paste0("taxonomy_", label, "_sequence_rownames.rds")))
}

message("[R] assignTaxonomy RDP16")
tax16 <- assignTaxonomy(
  seqs,
  rdp16,
  multithread = mt,
  tryRC = TRUE
)

format_tax(tax16, "rdp16_toGenus")

message("[R] assignTaxonomy RDP19")
tax19 <- assignTaxonomy(
  seqs,
  rdp19,
  multithread = mt,
  tryRC = TRUE
)

format_tax(tax19, "rdp19_toGenus")

message("[R] Done.")
message("[R] Main outputs:")
message("  ", file.path(results_dir, "asv_counts_by_run.tsv"))
message("  ", file.path(results_dir, "asv_counts_by_sample_accession.tsv"))
message("  ", file.path(results_dir, "asv_sequences.tsv"))
message("  ", file.path(results_dir, "taxonomy_rdp16_toGenus.tsv"))
message("  ", file.path(results_dir, "taxonomy_rdp19_toGenus.tsv"))
message("  ", file.path(results_dir, "read_tracking.tsv"))
RSCRIPT

############################################
# Run DADA2
############################################

echo "[dada2] Running R pipeline"

Rscript "${SCRIPT_DIR}/run_dada2_rdp_genus.R" \
  "${META_DIR}/trimmed_manifest.tsv" \
  "$RESULTS_DIR" \
  "$FILTER_DIR" \
  "${REF_DIR}/rdp_train_set_16.fa.gz" \
  "${REF_DIR}/rdp_19_toGenus_trainset.fa.gz" \
  "$THREADS" \
  "$TRUNC_F" \
  "$TRUNC_R" \
  "$MAXEE_F" \
  "$MAXEE_R" \
  "$MIN_OVERLAP"

echo "[job] Finished: $(date)"
echo "[job] Results:"
find "$RESULTS_DIR" -maxdepth 1 -type f -printf "%p\n" | sort
