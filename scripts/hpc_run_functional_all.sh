#!/usr/bin/env bash
set -Eeuo pipefail

PS4='+ $(date --iso-8601=seconds) ${BASH_SOURCE##*/}:${LINENO}: '
if [[ ${MUTT_TRACE:-0} == 1 ]]; then
  set -x
fi

log() {
  printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$*"
}

on_error() {
  local status=$?
  printf '[%s] ERROR status=%d script=%s line=%d command=%q\n' \
    "$(date --iso-8601=seconds)" "$status" "${BASH_SOURCE[1]##*/}" \
    "${BASH_LINENO[0]}" "$BASH_COMMAND" >&2
  exit "$status"
}
trap on_error ERR

readonly DEFAULT_REPOSITORY_URL="git@github.com:Silverman-Lab/mutt.git"
readonly DEFAULT_BRANCH="main"
readonly DEFAULT_CPUS="16"
readonly DEFAULT_MEMORY="64G"
readonly DEFAULT_IBD_TIME="2-00:00:00"
readonly DEFAULT_REST_TIME="2-00:00:00"
readonly DEFAULT_SETUP_TIME="02:00:00"
readonly DEFAULT_FINALIZE_TIME="04:00:00"
readonly DEFAULT_IBD_CONCURRENCY="2"
readonly DEFAULT_REST_CONCURRENCY="2"
readonly DEFAULT_MIN_FREE_GB="100"
readonly DEFAULT_MAX_LFS_BYTES="1900000000"

readonly -a IBD_STUDIES=(
  "2017_vandeputte_nature_flow"
  "2019_contijoch_elife_multispeciesqPCRshotgunandamplicon"
  "2019_vieirasilva_naturemicrobiology_pscibd"
  "2024_caenepeel_gastroenterology_IBD_CD_flow"
  "2024_sternes_frontmicrobiol_IBDppiqPCR"
)

readonly -a COVARIANCE_STUDIES=(
  "2016_stammler_microbiome_micehuman"
  "2017_liu_mbio_penilehivqPCR"
  "2017_vandeputte_nature_flow"
  "2019_contijoch_elife_multispeciesqPCRshotgunandamplicon"
  "2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein"
  "2019_morton_naturecommunications_songbird_oral"
  "2019_vieirasilva_naturemicrobiology_pscibd"
  "2020_barlow_naturecommunications_miceGI"
  "2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy"
  "2020_tettamantiboshier_msystems_vaginaltimeseries"
  "2020_zemb_microOpen_spike"
  "2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct"
  "2021_marotz_mSystems_oral_mouthwash"
  "2021_reese_cell_chimpanzee"
  "2021_vandeputte_naturecommunications_flow_timeseries"
  "2022_cvandevelde_ismecommunications_culturedflowhumanfecal"
  "2022_dreier_bmcmicrobiology_cheeseqpcr"
  "2022_jin_natureComm_technicalReplicates"
  "2022_krawczyk_microbiome_tickgeographicaldistributionqpcr"
  "2022_suriano_aps_micefecal"
  "2022_zaramela_msystems_synDNA"
  "2023_feng_imetawiley_chickensegment"
  "2023_kallastu_research_foodscience_food"
  "2023_maghini_naturebiotechnology_samplemesurement"
  "2023_pereira_nature_nervous"
  "2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota"
  "2024_caenepeel_gastroenterology_IBD_CD_flow"
  "2024_garciamartinez_bmcmicrobiology_ckdanddysbiosiswithserum"
  "2024_kruger_scientificreports_ddpcrhealthysubjects"
  "2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal"
  "2024_sternes_frontmicrobiol_IBDppiqPCR"
  "2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity"
  "2025_thiruppathy_microbiome_relicDNAflow"
)

usage() {
  cat <<'EOF'
Usage:
  SLURM_ACCOUNT=account ./hpc_run_functional_all.sh
  ./hpc_run_functional_all.sh --validate-only
  MUTT_PHASE=ibd|rest ./hpc_run_functional_all.sh --study-worker

The default command submits a setup job, one array task per IBD study, an IBD
finalizer, one array task per remaining study, and a final finalizer. Each stage
uses afterok dependencies. Required credentials are read from the user's Git or
SSH configuration; this script never stores credentials.

Useful environment overrides:
  MUTT_REPO_URL, MUTT_BRANCH, MUTT_WORK_ROOT, MUTT_IMAGE_URI
  MUTT_CPUS, MUTT_MEMORY, MUTT_IBD_TIME, MUTT_REST_TIME
  MUTT_SETUP_TIME, MUTT_FINALIZE_TIME
  MUTT_IBD_CONCURRENCY, MUTT_REST_CONCURRENCY
  MUTT_MIN_FREE_GB, MUTT_MAX_LFS_BYTES, SLURM_PARTITION
  APPTAINER_MODULE, MUTT_GIT_NAME, MUTT_GIT_EMAIL
  MUTT_FUNCTIONAL_MODE=use|rebuild|revalidate, MUTT_TRACE=0|1
EOF
}

contains_study() {
  local needle=$1
  shift
  local value
  for value in "$@"; do
    [[ "$value" == "$needle" ]] && return 0
  done
  return 1
}

validate_lists() {
  local -A seen=()
  local study
  for study in "${COVARIANCE_STUDIES[@]}"; do
    if [[ -n ${seen[$study]:-} ]]; then
      echo "Duplicate covariance study: $study" >&2
      return 1
    fi
    seen[$study]=1
  done
  for study in "${IBD_STUDIES[@]}"; do
    if ! contains_study "$study" "${COVARIANCE_STUDIES[@]}"; then
      echo "IBD study is absent from covariance cohort: $study" >&2
      return 1
    fi
  done
  [[ ${#IBD_STUDIES[@]} -eq 5 ]]
  [[ ${#COVARIANCE_STUDIES[@]} -eq 33 ]]
}

print_plan() {
  cat <<EOF
IBD studies: ${#IBD_STUDIES[@]}
Remaining studies: $((${#COVARIANCE_STUDIES[@]} - ${#IBD_STUDIES[@]}))
CPUs per task: ${MUTT_CPUS:-$DEFAULT_CPUS}
Memory per task: ${MUTT_MEMORY:-$DEFAULT_MEMORY}
IBD walltime: ${MUTT_IBD_TIME:-$DEFAULT_IBD_TIME}
Remaining walltime: ${MUTT_REST_TIME:-$DEFAULT_REST_TIME}
Walltime applies to each study array task.
IBD maximum concurrent tasks: ${MUTT_IBD_CONCURRENCY:-$DEFAULT_IBD_CONCURRENCY}
Remaining maximum concurrent tasks: ${MUTT_REST_CONCURRENCY:-$DEFAULT_REST_CONCURRENCY}
Setup walltime: ${MUTT_SETUP_TIME:-$DEFAULT_SETUP_TIME}
Finalizer walltime: ${MUTT_FINALIZE_TIME:-$DEFAULT_FINALIZE_TIME}
Minimum free scratch: ${MUTT_MIN_FREE_GB:-$DEFAULT_MIN_FREE_GB} GB
Functional mode: ${MUTT_FUNCTIONAL_MODE:-use}
EOF
}

submit_pipeline() {
  : "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT to the allocation that should be charged.}"
  command -v sbatch >/dev/null 2>&1 || {
    echo "sbatch was not found. Run this command from a Slurm login node." >&2
    exit 1
  }
  command -v git >/dev/null 2>&1 || { echo "git was not found." >&2; exit 1; }
  local concurrency
  for concurrency in \
      "${MUTT_IBD_CONCURRENCY:-$DEFAULT_IBD_CONCURRENCY}" \
      "${MUTT_REST_CONCURRENCY:-$DEFAULT_REST_CONCURRENCY}"; do
    [[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || {
      echo "Array concurrency limits must be positive integers." >&2
      exit 2
    }
  done

  local script repository_url branch code_commit run_id scratch work_root log_dir
  script=$(readlink -f "$0")
  repository_url=${MUTT_REPO_URL:-$DEFAULT_REPOSITORY_URL}
  branch=${MUTT_BRANCH:-$DEFAULT_BRANCH}
  code_commit=${MUTT_CODE_COMMIT:-}
  if [[ -z "$code_commit" ]]; then
    code_commit=$(git ls-remote "$repository_url" "refs/heads/$branch" | awk 'NR == 1 {print $1}')
  fi
  [[ "$code_commit" =~ ^[0-9a-f]{40}$ ]] || {
    echo "Could not resolve the remote $branch commit." >&2
    exit 1
  }
  run_id=${MUTT_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}
  scratch=${SCRATCH:-/scratch/$USER}
  work_root=${MUTT_WORK_ROOT:-$scratch/mutt-functional/$run_id}
  log_dir=$work_root/slurm_logs
  mkdir -p "$log_dir"

  local -a common=(
    --account="$SLURM_ACCOUNT"
  )
  if [[ -n ${SLURM_PARTITION:-} ]]; then
    common+=(--partition="$SLURM_PARTITION")
  fi

  local export_vars setup_job ibd_job ibd_finalize_job rest_job rest_finalize_job
  local ibd_last rest_last
  export_vars="ALL,MUTT_RUN_ID=$run_id,MUTT_WORK_ROOT=$work_root,MUTT_CODE_COMMIT=$code_commit,MUTT_REPO_URL=$repository_url,MUTT_BRANCH=$branch"
  ibd_last=$((${#IBD_STUDIES[@]} - 1))
  rest_last=$((${#COVARIANCE_STUDIES[@]} - ${#IBD_STUDIES[@]} - 1))

  setup_job=$(sbatch --parsable "${common[@]}" \
    --job-name=mutt_setup \
    --cpus-per-task=2 \
    --mem=8G \
    --time="${MUTT_SETUP_TIME:-$DEFAULT_SETUP_TIME}" \
    --output="$log_dir/%x_%j.out" \
    --error="$log_dir/%x_%j.err" \
    --export="$export_vars" \
    "$script" --setup-worker)
  setup_job=${setup_job%%;*}
  ibd_job=$(sbatch --parsable "${common[@]}" \
    --dependency="afterok:$setup_job" \
    --job-name=mutt_ibd \
    --array="0-${ibd_last}%${MUTT_IBD_CONCURRENCY:-$DEFAULT_IBD_CONCURRENCY}" \
    --cpus-per-task="${MUTT_CPUS:-$DEFAULT_CPUS}" \
    --mem="${MUTT_MEMORY:-$DEFAULT_MEMORY}" \
    --time="${MUTT_IBD_TIME:-$DEFAULT_IBD_TIME}" \
    --output="$log_dir/%x_%A_%a.out" \
    --error="$log_dir/%x_%A_%a.err" \
    --export="$export_vars,MUTT_PHASE=ibd" \
    "$script" --study-worker)
  ibd_job=${ibd_job%%;*}
  ibd_finalize_job=$(sbatch --parsable "${common[@]}" \
    --dependency="afterok:$ibd_job" \
    --job-name=mutt_ibd_finalize \
    --cpus-per-task=2 \
    --mem=8G \
    --time="${MUTT_FINALIZE_TIME:-$DEFAULT_FINALIZE_TIME}" \
    --output="$log_dir/%x_%j.out" \
    --error="$log_dir/%x_%j.err" \
    --export="$export_vars,MUTT_PHASE=ibd" \
    "$script" --finalize-worker)
  ibd_finalize_job=${ibd_finalize_job%%;*}
  rest_job=$(sbatch --parsable "${common[@]}" \
    --dependency="afterok:$ibd_finalize_job" \
    --job-name=mutt_rest \
    --array="0-${rest_last}%${MUTT_REST_CONCURRENCY:-$DEFAULT_REST_CONCURRENCY}" \
    --cpus-per-task="${MUTT_CPUS:-$DEFAULT_CPUS}" \
    --mem="${MUTT_MEMORY:-$DEFAULT_MEMORY}" \
    --time="${MUTT_REST_TIME:-$DEFAULT_REST_TIME}" \
    --output="$log_dir/%x_%A_%a.out" \
    --error="$log_dir/%x_%A_%a.err" \
    --export="$export_vars,MUTT_PHASE=rest" \
    "$script" --study-worker)
  rest_job=${rest_job%%;*}
  rest_finalize_job=$(sbatch --parsable "${common[@]}" \
    --dependency="afterok:$rest_job" \
    --job-name=mutt_rest_finalize \
    --cpus-per-task=2 \
    --mem=8G \
    --time="${MUTT_FINALIZE_TIME:-$DEFAULT_FINALIZE_TIME}" \
    --output="$log_dir/%x_%j.out" \
    --error="$log_dir/%x_%j.err" \
    --export="$export_vars,MUTT_PHASE=rest" \
    "$script" --finalize-worker)
  rest_finalize_job=${rest_finalize_job%%;*}

  print_plan
  cat <<EOF
Run ID: $run_id
Pinned package commit: $code_commit
Work root: $work_root
Setup job: $setup_job
IBD array: $ibd_job (afterok:$setup_job)
IBD finalizer: $ibd_finalize_job (afterok:$ibd_job)
Remaining-studies array: $rest_job (afterok:$ibd_finalize_job)
Remaining-studies finalizer: $rest_finalize_job (afterok:$rest_job)
Monitor: squeue -j $setup_job,$ibd_job,$ibd_finalize_job,$rest_job,$rest_finalize_job
After completion: sacct -j $setup_job,$ibd_job,$ibd_finalize_job,$rest_job,$rest_finalize_job --format=JobID,State,Elapsed,MaxRSS,ExitCode
EOF
}

load_apptainer() {
  if ! command -v apptainer >/dev/null 2>&1 && [[ -n ${APPTAINER_MODULE:-} ]] &&
      command -v module >/dev/null 2>&1; then
    module load "$APPTAINER_MODULE"
  fi
  command -v apptainer >/dev/null 2>&1 || {
    echo "apptainer was not found. Set APPTAINER_MODULE if the cluster provides a module." >&2
    exit 1
  }
}

check_free_space() {
  local path=$1 minimum_gb=$2 available_kb required_kb
  available_kb=$(df -Pk "$path" | awk 'NR == 2 {print $4}')
  required_kb=$((minimum_gb * 1024 * 1024))
  if (( available_kb < required_kb )); then
    echo "Insufficient free space under $path: require at least ${minimum_gb} GB." >&2
    exit 1
  fi
}

prepare_repository() {
  local repository_url=$1 branch=$2 repository=$3
  if [[ ! -d "$repository/.git" ]]; then
    GIT_LFS_SKIP_SMUDGE=1 git clone --branch "$branch" "$repository_url" "$repository"
  fi
  cd "$repository"
  [[ $(git remote get-url origin) == "$repository_url" ]] || {
    echo "Existing checkout has an unexpected origin: $repository" >&2
    exit 1
  }
  if [[ -n $(git status --porcelain --untracked-files=all) ]]; then
    echo "Tracked or untracked changes already exist in $repository; refusing to overwrite them." >&2
    git status --short >&2
    exit 1
  fi
  git lfs install --local
  git fetch origin "$branch"
  git checkout "$branch"
  git pull --ff-only origin "$branch"
  git lfs pull
}

prepare_image() {
  local work_root=$1 code_commit=$2 repository=$3
  local image_dir image temporary uri
  image_dir=$work_root/images
  image=$image_dir/mutt-${code_commit}.sif
  mkdir -p "$image_dir"
  if [[ ! -f "$image" ]]; then
    uri=${MUTT_IMAGE_URI:-docker://ghcr.io/silverman-lab/mutt:sha-$code_commit}
    temporary=${image}.partial.${SLURM_JOB_ID:-$$}
    log "Pulling Apptainer image: $uri" >&2
    if ! apptainer pull "$temporary" "$uri" >&2; then
      echo "Apptainer could not pull $uri. For private GHCR images, run 'apptainer registry login docker://ghcr.io' on the HPC login node or set MUTT_IMAGE_URI to a readable image." >&2
      return 1
    fi
    [[ -s "$temporary" ]] || {
      echo "Apptainer reported success but did not create $temporary." >&2
      return 1
    }
    mv "$temporary" "$image"
  fi
  echo "$image"
}

write_checkpoint() {
  local marker=$1 status=$2 temporary
  mkdir -p "$(dirname "$marker")"
  temporary=${marker}.partial.${SLURM_JOB_ID:-$$}
  {
    echo "status=$status"
    echo "run_id=$MUTT_RUN_ID"
    echo "code_commit=$MUTT_CODE_COMMIT"
    echo "job_id=${SLURM_JOB_ID:-NA}"
    echo "completed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "$temporary"
  mv "$temporary" "$marker"
  log "Checkpoint written: $marker"
}

install_source_package() {
  local image=$1 repository=$2 library=$3
  mkdir -p "$library"
  apptainer exec --cleanenv \
    --bind "$repository:/work" \
    --bind "$library:/mutt-lib" \
    --pwd /work \
    "$image" \
    R CMD INSTALL --library=/mutt-lib /work
  apptainer exec --cleanenv \
    --bind "$repository:/work" \
    --bind "$library:/mutt-lib" \
    --pwd /work \
    --env R_LIBS_USER=/mutt-lib \
    "$image" \
    Rscript --vanilla -e 'library(mutt); stopifnot(identical(getNamespaceExports("mutt"), "mutt")); cat("MUTT source installation passed\n")'
}

phase_studies() {
  local phase=$1 study
  if [[ "$phase" == "ibd" ]]; then
    printf '%s\n' "${IBD_STUDIES[@]}"
    return
  fi
  for study in "${COVARIANCE_STUDIES[@]}"; do
    if ! contains_study "$study" "${IBD_STUDIES[@]}"; then
      printf '%s\n' "$study"
    fi
  done
}

phase_study_at_index() {
  local phase=$1 index=$2
  mapfile -t studies < <(phase_studies "$phase")
  [[ "$index" =~ ^[0-9]+$ ]] && (( index < ${#studies[@]} )) || {
    echo "Invalid array index $index for phase $phase." >&2
    exit 2
  }
  printf '%s\n' "${studies[$index]}"
}

validate_registry() {
  local repository=$1 study
  while IFS= read -r study; do
    if ! awk -F, -v target="$study" 'NR > 1 && $1 == target {found=1} END {exit !found}' \
        "$repository/inst/parsers/index.csv"; then
      echo "Study is not registered in MUTT: $study" >&2
      exit 1
    fi
    [[ -d "$repository/studies/$study" ]] || {
      echo "Study directory is missing after Git LFS pull: $study" >&2
      exit 1
    }
  done < <(phase_studies "$MUTT_PHASE")
}

run_study() {
  local image=$1 repository=$2 library=$3 phase_output=$4 study=$5 mode=$6 cpus=$7
  echo "[$(date --iso-8601=seconds)] Starting $study"
  apptainer exec --cleanenv \
    --bind "$repository:/work" \
    --bind "$library:/mutt-lib" \
    --bind "$phase_output:/phase-output" \
    --pwd /work \
    --env R_LIBS_USER=/mutt-lib \
    --env MUTT_DATA_DIR=/work/studies \
    --env MUTT_FUNCTIONAL_PROCESSES="$cpus" \
    "$image" \
    Rscript --vanilla /work/scripts/run_functional_study.R \
      "$study" /work /phase-output "$mode"
  echo "[$(date --iso-8601=seconds)] Completed $study"
}

# Retained for recovery of runs created before functional_data schema version 1.
archive_study_cache_legacy() {
  local repository=$1 study=$2 max_bytes=$3 archive bytes pointer
  archive="$repository/studies/$study/functional.zip"
  [[ -d "$repository/studies/$study/functional" ]] || {
    echo "Functional cache directory was not created for $study." >&2
    exit 1
  }
  (
    cd "$repository"
    ./zip-push-gitlfs.sh --apply --keep-source --replace "studies/$study/functional"
    git add -- "studies/$study/functional.zip"
  )
  bytes=$(stat -c '%s' "$archive")
  if (( bytes > max_bytes )); then
    echo "Archive exceeds configured Git LFS per-file limit ($max_bytes bytes): $archive" >&2
    exit 1
  fi
  pointer=$(git -C "$repository" show ":studies/$study/functional.zip" | sed -n '1p')
  [[ "$pointer" == "version https://git-lfs.github.com/spec/v1" ]] || {
    echo "Archive was not staged as a Git LFS pointer: $archive" >&2
    exit 1
  }
}

archive_study_cache() {
    local repository=$1 study=$2 max_bytes=$3 publication file relative pointer
    publication="$repository/studies/$study/functional_data"

    [[ -d "$repository/studies/$study/functional" ]] || {
        echo "Functional execution cache was not created for $study." >&2
        exit 1
    }

    log "Building adaptive functional bundle for $study."
    apptainer exec --cleanenv \
        --bind "$repository:/work" \
        --bind "$library:/mutt-lib" \
        --pwd /work \
        --env R_LIBS_USER=/mutt-lib \
        "$image" \
        Rscript --vanilla /work/scripts/build_functional_bundle.R \
        --study-dir "/work/studies/$study" \
        --max-asset-bytes "$max_bytes"

    [[ -f "$publication/bundle_manifest.json" ]] || {
        echo "Functional bundle manifest was not created for $study." >&2
        exit 1
    }
    [[ -f "$publication/functional-core.zip" ]] || {
        echo "Functional core archive was not created for $study." >&2
        exit 1
    }

    while IFS= read -r -d '' file; do
        if (( $(stat -c '%s' "$file") > max_bytes )); then
            echo "Published functional asset exceeds limit ($max_bytes bytes): $file" >&2
            exit 1
        fi
    done < <(find "$publication" -type f -print0)

    git -C "$repository" add -- "studies/$study/functional_data"
    while IFS= read -r -d '' file; do
        [[ "$file" == *.json ]] && continue
        relative=${file#"$repository/"}
        pointer=$(git -C "$repository" show ":$relative" | sed -n '1p')
        [[ "$pointer" == "version https://git-lfs.github.com/spec/v1" ]] || {
            echo "Published functional asset is not staged through Git LFS: $relative" >&2
            exit 1
        }
    done < <(find "$publication" -type f -print0)
}

commit_and_push_phase() {
  local repository=$1 phase=$2 run_id=$3 branch=$4 output_dir=$5
  cd "$repository"
  git add -- "$output_dir"
  git diff --cached --check
  if git diff --cached --quiet; then
    echo "No phase outputs were available to commit." >&2
    exit 1
  fi
  if [[ -n ${MUTT_GIT_NAME:-} ]]; then git config user.name "$MUTT_GIT_NAME"; fi
  if [[ -n ${MUTT_GIT_EMAIL:-} ]]; then git config user.email "$MUTT_GIT_EMAIL"; fi
  git config user.name >/dev/null
  git config user.email >/dev/null
  git commit -m "Add MUTT functional outputs: $phase ($run_id)"
  git fetch origin "$branch"
  git rebase "origin/$branch"
  git push origin "HEAD:$branch"
}

validate_worker_environment() {
  : "${MUTT_RUN_ID:?MUTT_RUN_ID is required.}"
  : "${MUTT_WORK_ROOT:?MUTT_WORK_ROOT is required.}"
  : "${MUTT_CODE_COMMIT:?MUTT_CODE_COMMIT is required.}"
}

validate_phase() {
  : "${MUTT_PHASE:?MUTT_PHASE must be ibd or rest.}"
  [[ "$MUTT_PHASE" == "ibd" || "$MUTT_PHASE" == "rest" ]] || {
    echo "MUTT_PHASE must be ibd or rest." >&2
    exit 2
  }
}

print_worker_settings() {
  local role=$1 repository_url=$2
  echo "Job ID: ${SLURM_JOB_ID:-NA}"
  echo "Host: $(hostname)"
  echo "Role: $role"
  echo "Phase: ${MUTT_PHASE:-NA}"
  echo "Array task: ${SLURM_ARRAY_TASK_ID:-NA}"
  echo "Run ID: $MUTT_RUN_ID"
  echo "Work root: $MUTT_WORK_ROOT"
  echo "Repository: $repository_url"
  echo "Pinned code commit: $MUTT_CODE_COMMIT"
  echo "Started: $(date --iso-8601=seconds)"
}

run_setup_worker() {
  validate_worker_environment

  local repository_url branch repository min_free image library current_head checkpoint
  repository_url=${MUTT_REPO_URL:-$DEFAULT_REPOSITORY_URL}
  branch=${MUTT_BRANCH:-$DEFAULT_BRANCH}
  repository=$MUTT_WORK_ROOT/repository
  min_free=${MUTT_MIN_FREE_GB:-$DEFAULT_MIN_FREE_GB}
  checkpoint=$MUTT_WORK_ROOT/checkpoints/setup.complete

  mkdir -p "$MUTT_WORK_ROOT/slurm_logs"
  if [[ -f "$checkpoint" ]]; then
    log "Setup checkpoint already exists; validating reusable setup."
    image=$MUTT_WORK_ROOT/images/mutt-${MUTT_CODE_COMMIT}.sif
    library=$MUTT_WORK_ROOT/r-library/$MUTT_CODE_COMMIT
    [[ -d "$repository/.git" && -s "$image" && -d "$library/mutt" ]] || {
      echo "Setup checkpoint exists, but one or more setup artifacts are missing." >&2
      exit 1
    }
    log "Setup artifacts are complete; skipping setup."
    return
  fi
  check_free_space "$MUTT_WORK_ROOT" "$min_free"
  load_apptainer
  command -v git >/dev/null 2>&1
  git lfs version
  print_worker_settings setup "$repository_url"
  command -v quota >/dev/null 2>&1 && quota -s || true

  log "Preparing repository checkout."
  prepare_repository "$repository_url" "$branch" "$repository"
  current_head=$(git -C "$repository" rev-parse HEAD)
  [[ "$current_head" == "$MUTT_CODE_COMMIT" ]] || {
    echo "Remote $branch changed after submission; resubmit to pin the new commit." >&2
    exit 1
  }
  log "Preparing container image."
  if ! image=$(prepare_image "$MUTT_WORK_ROOT" "$MUTT_CODE_COMMIT" "$repository"); then
    echo "Container image preparation failed." >&2
    exit 1
  fi
  library=$MUTT_WORK_ROOT/r-library/$MUTT_CODE_COMMIT
  echo "Image: $image"
  log "Installing and smoke-testing the MUTT source package."
  install_source_package "$image" "$repository" "$library"
  log "Validating the IBD and remaining-study registries."
  MUTT_PHASE=ibd validate_registry "$repository"
  MUTT_PHASE=rest validate_registry "$repository"
  write_checkpoint "$checkpoint" setup_complete
  echo "Setup completed: $(date --iso-8601=seconds)"
}

run_study_worker() {
  validate_worker_environment
  validate_phase
  : "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required.}"

  local repository_url repository mode cpus image library phase_output study checkpoint
  repository_url=${MUTT_REPO_URL:-$DEFAULT_REPOSITORY_URL}
  repository=$MUTT_WORK_ROOT/repository
  mode=${MUTT_FUNCTIONAL_MODE:-use}
  cpus=${SLURM_CPUS_PER_TASK:-${MUTT_CPUS:-$DEFAULT_CPUS}}
  case "$mode" in
    use|rebuild|revalidate)
      ;;
    *)
      echo "MUTT_FUNCTIONAL_MODE must be use, rebuild, or revalidate." >&2
      exit 2
      ;;
  esac
  [[ -f "$MUTT_WORK_ROOT/checkpoints/setup.complete" ]] || {
    echo "Setup completion marker is missing." >&2
    exit 1
  }
  load_apptainer
  print_worker_settings study "$repository_url"
  echo "CPUs: $cpus"
  echo "Functional mode: $mode"

  image=$MUTT_WORK_ROOT/images/mutt-${MUTT_CODE_COMMIT}.sif
  library=$MUTT_WORK_ROOT/r-library/$MUTT_CODE_COMMIT
  [[ -f "$image" && -d "$library/mutt" ]] || {
    echo "Prepared image or R package library is missing." >&2
    exit 1
  }
  echo "Image: $image"
  phase_output="$MUTT_WORK_ROOT/phase-output/$MUTT_PHASE"
  mkdir -p "$phase_output"
  study=$(phase_study_at_index "$MUTT_PHASE" "$SLURM_ARRAY_TASK_ID")
  checkpoint=$MUTT_WORK_ROOT/checkpoints/$MUTT_PHASE/$study.complete
  if [[ -f "$checkpoint" ]]; then
    log "Study checkpoint already exists; skipping $study."
    return
  fi
  log "Running study $study."
  run_study "$image" "$repository" "$library" "$phase_output" "$study" "$mode" "$cpus"
  write_checkpoint "$checkpoint" study_complete
}

run_finalize_worker() {
  validate_worker_environment
  validate_phase

  local repository_url branch repository max_lfs output_dir phase_output study checkpoint
  repository_url=${MUTT_REPO_URL:-$DEFAULT_REPOSITORY_URL}
  branch=${MUTT_BRANCH:-$DEFAULT_BRANCH}
  repository=$MUTT_WORK_ROOT/repository
  max_lfs=${MUTT_MAX_LFS_BYTES:-$DEFAULT_MAX_LFS_BYTES}
  output_dir="analysis/hpc_functional/$MUTT_RUN_ID/$MUTT_PHASE"
  phase_output="$MUTT_WORK_ROOT/phase-output/$MUTT_PHASE"
  checkpoint=$MUTT_WORK_ROOT/checkpoints/$MUTT_PHASE/finalize.complete

  print_worker_settings finalize "$repository_url"
  if [[ -f "$checkpoint" ]]; then
    log "Finalization checkpoint already exists; skipping phase $MUTT_PHASE."
    return
  fi
  command -v git >/dev/null 2>&1
  git lfs version
  command -v zip >/dev/null 2>&1
  command -v unzip >/dev/null 2>&1

  while IFS= read -r study; do
    [[ -f "$MUTT_WORK_ROOT/checkpoints/$MUTT_PHASE/$study.complete" ]] || {
      echo "Study checkpoint is missing for $study; refusing to finalize." >&2
      exit 1
    }
        log "Publishing API-facing functional data for $study."
        archive_study_cache "$repository" "$study" "$max_lfs"
  done < <(phase_studies "$MUTT_PHASE")

  {
    echo "phase=$MUTT_PHASE"
    echo "run_id=$MUTT_RUN_ID"
    echo "code_commit=$MUTT_CODE_COMMIT"
    echo "completed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "$phase_output/phase_complete.txt"
  mkdir -p "$repository/$output_dir"
  cp -a "$phase_output/." "$repository/$output_dir/"
  log "Committing and pushing phase $MUTT_PHASE."
  commit_and_push_phase "$repository" "$MUTT_PHASE" "$MUTT_RUN_ID" "$branch" "$output_dir"
  write_checkpoint "$checkpoint" finalize_complete
  echo "Finished and pushed: $(date --iso-8601=seconds)"
}

validate_lists
case ${1:-} in
  --validate-only)
    print_plan
    printf 'IBD cohort:\n%s\n' "$(printf '  %s\n' "${IBD_STUDIES[@]}")"
    echo "Validation passed. No jobs were submitted."
    ;;
  --setup-worker)
    run_setup_worker
    ;;
  --study-worker)
    run_study_worker
    ;;
  --finalize-worker)
    run_finalize_worker
    ;;
  -h|--help)
    usage
    ;;
  "")
    submit_pipeline
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
