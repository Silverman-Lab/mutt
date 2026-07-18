#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-dada2_rdp}"
PROJECT_DIR="${PROJECT_DIR:-$PWD/PRJEB71738_dada2_rdp}"
REF_DIR="${PROJECT_DIR}/refs"

mkdir -p "$REF_DIR"

echo "[setup] Project dir: $PROJECT_DIR"
echo "[setup] Reference dir: $REF_DIR"
echo "[setup] Conda env: $ENV_NAME"

if command -v mamba >/dev/null 2>&1; then
  CONDA_FRONTEND="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_FRONTEND="conda"
else
  echo "ERROR: neither conda nor mamba is available on PATH." >&2
  echo "Load your HPC conda/miniforge module first, then rerun." >&2
  exit 1
fi

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[setup] Conda env already exists: $ENV_NAME"
else
  echo "[setup] Creating conda env..."
  "$CONDA_FRONTEND" create -y -n "$ENV_NAME" \
    -c conda-forge -c bioconda \
    r-base \
    bioconductor-dada2 \
    cutadapt \
    curl \
    wget \
    python \
    pigz \
    coreutils
fi

download_if_missing() {
  local url="$1"
  local out="$2"

  if [[ -s "$out" ]]; then
    echo "[setup] Already exists: $out"
  else
    echo "[setup] Downloading: $out"
    curl -L --fail --retry 5 --retry-delay 10 -o "$out" "$url"
  fi
}

# RDP trainset 16, genus-level.
download_if_missing \
  "https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz?download=1" \
  "${REF_DIR}/rdp_train_set_16.fa.gz"

# RDP trainset 19, genus-level.
download_if_missing \
  "https://zenodo.org/records/14168771/files/rdp_19_toGenus_trainset.fa.gz?download=1" \
  "${REF_DIR}/rdp_19_toGenus_trainset.fa.gz"

echo "[setup] Checking files:"
ls -lh "$REF_DIR"

echo "[setup] Done."
echo
echo "Next:"
echo "  sbatch 01_run_PRJEB71738_dada2_rdp.sbatch"
