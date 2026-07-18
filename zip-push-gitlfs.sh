#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./zip-push-gitlfs.sh [--apply] [--keep-source] [--replace] STUDY_PATH [...]

Create one ZIP archive per explicitly listed study-data file, or archive a
study's `functional/` cache directory. The default is a dry run. Pass --apply
to create validated archives. Files are removed after validation unless
--keep-source is supplied; directories require --keep-source. Pass --replace
to atomically refresh an existing archive. This script never changes remotes,
stages, commits, or pushes Git changes.
EOF
}

apply=false
keep_source=false
replace=false
files=()
for argument in "$@"; do
  case "$argument" in
    --apply) apply=true ;;
    --keep-source) keep_source=true ;;
    --replace) replace=true ;;
    -h|--help) usage; exit 0 ;;
    --*) echo "Unknown option: $argument" >&2; usage >&2; exit 2 ;;
    *) files+=("$argument") ;;
  esac
done

if ! command -v git >/dev/null 2>&1 || ! command -v git-lfs >/dev/null 2>&1; then
  echo "Git and Git LFS are required." >&2
  exit 1
fi
if ! command -v zip >/dev/null 2>&1 || ! command -v unzip >/dev/null 2>&1 ||
    ! command -v diff >/dev/null 2>&1; then
  echo "zip, unzip, and diff are required." >&2
  exit 1
fi
if [[ ${#files[@]} -eq 0 ]]; then
  usage >&2
  exit 2
fi

repository=$(git rev-parse --show-toplevel)
cd "$repository"

for path in "${files[@]}"; do
  path=${path#./}
  path=${path%/}
  if [[ "$path" != studies/* ]]; then
    echo "Refusing to archive a path outside studies/: $path" >&2
    exit 1
  fi
  if [[ ! -e "$path" || -L "$path" ]]; then
    echo "Study-data path does not exist or is a symbolic link: $path" >&2
    exit 1
  fi
  if [[ -d "$path" ]]; then
    if [[ $(basename "$path") != "functional" ]]; then
      echo "Only study functional/ directories may be archived: $path" >&2
      exit 1
    fi
    if [[ "$keep_source" != true ]]; then
      echo "Directory archives require --keep-source: $path" >&2
      exit 1
    fi
    if [[ -n $(find "$path" -type l -print -quit) ]]; then
      echo "Refusing to archive a directory containing symbolic links: $path" >&2
      exit 1
    fi
  else
    case "$path" in
      *.zip|*.gz|*.R|*.Rmd|*.sh|*.sbatch|*.py|*.cpp|*.md)
        echo "Refusing to archive source, documentation, or an existing archive: $path" >&2
        exit 1
        ;;
    esac
  fi

  archive="${path}.zip"
  if [[ -e "$archive" && "$replace" != true ]]; then
    echo "Refusing to overwrite existing archive: $archive" >&2
    exit 1
  fi
  filter=$(git check-attr filter -- "$archive" | sed 's/^.*: filter: //')
  if [[ "$filter" != "lfs" ]]; then
    echo "Archive is not covered by Git LFS attributes: $archive" >&2
    exit 1
  fi

  if [[ "$apply" == false ]]; then
    action="create"
    [[ -e "$archive" ]] && action="replace"
    echo "Would $action: $path -> $archive"
    continue
  fi

  temporary="${archive}.tmp.$$.zip"
  extracted=$(mktemp -d "${TMPDIR:-/tmp}/mutt-zip-verify-XXXXXX")
  trap 'rm -f -- "$temporary"; rm -rf -- "$extracted"' EXIT
  if [[ -d "$path" ]]; then
    parent=$(dirname "$path")
    name=$(basename "$path")
    absolute_temporary="$repository/$temporary"
    (cd "$parent" && zip -q -X -r "$absolute_temporary" "$name")
    unzip -tq "$temporary" >/dev/null
    unzip -q "$temporary" -d "$extracted"
    diff -qr --no-dereference "$path" "$extracted/$name" >/dev/null
  else
    zip -q -X -j "$temporary" "$path"
    unzip -tq "$temporary" >/dev/null
    unzip -p "$temporary" "$(basename "$path")" | cmp - "$path"
  fi
  if [[ "$replace" == true ]]; then
    mv -f -- "$temporary" "$archive"
  else
    mv -- "$temporary" "$archive"
  fi
  if [[ "$keep_source" != true ]]; then
    rm -- "$path"
  fi
  rm -rf -- "$extracted"
  trap - EXIT
  echo "Archived: $path -> $archive"
done

if [[ "$apply" == false ]]; then
  echo "Dry run only. Re-run with --apply after reviewing this list."
else
  echo "Archives validated. Review with: git status --short && git lfs status"
fi
