#!/usr/bin/env bash
set -euo pipefail

readonly DEFAULT_REPOSITORY_URL="git@github.com:Silverman-Lab/mutt.git"
readonly DEFAULT_BRANCH="main"
readonly DEFAULT_CPUS="16"
readonly DEFAULT_MEMORY="64G"
readonly DEFAULT_IBD_TIME="2-00:00:00"
readonly DEFAULT_REST_TIME="7-00:00:00"
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
  MUTT_PHASE=ibd|rest ./hpc_run_functional_all.sh --worker

The default command submits two Slurm jobs. The rest job has an afterok
dependency on the IBD job. Required credentials are read from the user's Git or
SSH configuration; this script never stores credentials.

Useful environment overrides:
  MUTT_REPO_URL, MUTT_BRANCH, MUTT_WORK_ROOT, MUTT_IMAGE_URI
  MUTT_CPUS, MUTT_MEMORY, MUTT_IBD_TIME, MUTT_REST_TIME
  MUTT_MIN_FREE_GB, MUTT_MAX_LFS_BYTES, SLURM_PARTITION
  APPTAINER_MODULE, MUTT_GIT_NAME, MUTT_GIT_EMAIL
  MUTT_FUNCTIONAL_MODE=use|rebuild
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
    --cpus-per-task="${MUTT_CPUS:-$DEFAULT_CPUS}"
    --mem="${MUTT_MEMORY:-$DEFAULT_MEMORY}"
    --output="$log_dir/%x_%j.out"
    --error="$log_dir/%x_%j.err"
  )
  if [[ -n ${SLURM_PARTITION:-} ]]; then
    common+=(--partition="$SLURM_PARTITION")
  fi

  local ibd_job rest_job
  ibd_job=$(sbatch --parsable "${common[@]}" \
    --job-name=mutt_ibd \
    --time="${MUTT_IBD_TIME:-$DEFAULT_IBD_TIME}" \
    --export="ALL,MUTT_PHASE=ibd,MUTT_RUN_ID=$run_id,MUTT_WORK_ROOT=$work_root,MUTT_CODE_COMMIT=$code_commit,MUTT_REPO_URL=$repository_url,MUTT_BRANCH=$branch" \
    "$script" --worker)
  ibd_job=${ibd_job%%;*}
  rest_job=$(sbatch --parsable "${common[@]}" \
    --dependency="afterok:$ibd_job" \
    --job-name=mutt_rest \
    --time="${MUTT_REST_TIME:-$DEFAULT_REST_TIME}" \
    --export="ALL,MUTT_PHASE=rest,MUTT_RUN_ID=$run_id,MUTT_WORK_ROOT=$work_root,MUTT_CODE_COMMIT=$code_commit,MUTT_REPO_URL=$repository_url,MUTT_BRANCH=$branch" \
    "$script" --worker)
  rest_job=${rest_job%%;*}

  print_plan
  cat <<EOF
Run ID: $run_id
Pinned package commit: $code_commit
Work root: $work_root
IBD job: $ibd_job
Remaining-studies job: $rest_job (afterok:$ibd_job)
Monitor: squeue -j $ibd_job,$rest_job
After completion: sacct -j $ibd_job,$rest_job --format=JobID,State,Elapsed,MaxRSS,ExitCode
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
    temporary=${image}.partial
    echo "Pulling Apptainer image: $uri" >&2
    apptainer pull "$temporary" "$uri" >&2
    mv "$temporary" "$image"
  fi
  echo "$image"
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

archive_study_cache() {
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

run_worker() {
  : "${MUTT_PHASE:?MUTT_PHASE must be ibd or rest.}"
  [[ "$MUTT_PHASE" == "ibd" || "$MUTT_PHASE" == "rest" ]] || {
    echo "MUTT_PHASE must be ibd or rest." >&2
    exit 2
  }
  : "${MUTT_RUN_ID:?MUTT_RUN_ID is required.}"
  : "${MUTT_WORK_ROOT:?MUTT_WORK_ROOT is required.}"
  : "${MUTT_CODE_COMMIT:?MUTT_CODE_COMMIT is required.}"

  local repository_url branch repository mode cpus min_free max_lfs image library
  local output_dir phase_output study current_head
  repository_url=${MUTT_REPO_URL:-$DEFAULT_REPOSITORY_URL}
  branch=${MUTT_BRANCH:-$DEFAULT_BRANCH}
  repository=$MUTT_WORK_ROOT/repository
  mode=${MUTT_FUNCTIONAL_MODE:-use}
  cpus=${SLURM_CPUS_PER_TASK:-${MUTT_CPUS:-$DEFAULT_CPUS}}
  min_free=${MUTT_MIN_FREE_GB:-$DEFAULT_MIN_FREE_GB}
  max_lfs=${MUTT_MAX_LFS_BYTES:-$DEFAULT_MAX_LFS_BYTES}
  [[ "$mode" == "use" || "$mode" == "rebuild" ]] || {
    echo "MUTT_FUNCTIONAL_MODE must be use or rebuild." >&2
    exit 2
  }

  mkdir -p "$MUTT_WORK_ROOT/slurm_logs"
  exec 9>"$MUTT_WORK_ROOT/pipeline.lock"
  flock -n 9 || { echo "Another MUTT phase is using $MUTT_WORK_ROOT." >&2; exit 1; }
  check_free_space "$MUTT_WORK_ROOT" "$min_free"
  load_apptainer
  command -v git >/dev/null 2>&1
  git lfs version
  command -v flock >/dev/null 2>&1
  command -v zip >/dev/null 2>&1
  command -v unzip >/dev/null 2>&1

  echo "Job ID: ${SLURM_JOB_ID:-NA}"
  echo "Host: $(hostname)"
  echo "Phase: $MUTT_PHASE"
  echo "Run ID: $MUTT_RUN_ID"
  echo "Work root: $MUTT_WORK_ROOT"
  echo "Repository: $repository_url"
  echo "Pinned code commit: $MUTT_CODE_COMMIT"
  echo "CPUs: $cpus"
  echo "Functional mode: $mode"
  echo "Started: $(date --iso-8601=seconds)"
  command -v quota >/dev/null 2>&1 && quota -s || true

  prepare_repository "$repository_url" "$branch" "$repository"
  current_head=$(git -C "$repository" rev-parse HEAD)
  if [[ "$MUTT_PHASE" == "ibd" && "$current_head" != "$MUTT_CODE_COMMIT" ]]; then
    echo "Remote main changed after submission; resubmit to pin the new code commit." >&2
    exit 1
  fi
  git -C "$repository" merge-base --is-ancestor "$MUTT_CODE_COMMIT" "$current_head" || {
    echo "Pinned package commit is not an ancestor of the phase checkout." >&2
    exit 1
  }

  image=$(prepare_image "$MUTT_WORK_ROOT" "$MUTT_CODE_COMMIT" "$repository")
  library=$MUTT_WORK_ROOT/r-library/$MUTT_CODE_COMMIT
  echo "Image: $image"
  install_source_package "$image" "$repository" "$library"
  validate_registry "$repository"

  output_dir="analysis/hpc_functional/$MUTT_RUN_ID/$MUTT_PHASE"
  phase_output="$MUTT_WORK_ROOT/phase-output/$MUTT_PHASE"
  mkdir -p "$phase_output"
  while IFS= read -r study; do
    run_study "$image" "$repository" "$library" "$phase_output" "$study" "$mode" "$cpus"
  done < <(phase_studies "$MUTT_PHASE")

  while IFS= read -r study; do
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
  commit_and_push_phase "$repository" "$MUTT_PHASE" "$MUTT_RUN_ID" "$branch" "$output_dir"
  echo "Finished and pushed: $(date --iso-8601=seconds)"
}

validate_lists
case ${1:-} in
  --validate-only)
    print_plan
    printf 'IBD cohort:\n%s\n' "$(printf '  %s\n' "${IBD_STUDIES[@]}")"
    echo "Validation passed. No jobs were submitted."
    ;;
  --worker)
    run_worker
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
