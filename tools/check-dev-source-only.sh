#!/usr/bin/env bash

set -euo pipefail
shopt -s nocasematch

readonly max_source_bytes=$((512 * 1024))
status=0

while IFS= read -r -d '' path; do
  case "$path" in
    *.png|*.jpg|*.jpeg|*.gif|*.webp|*.pdf|*.html|*.rds|*.rdata|*.csv)
      printf 'Generated artifact is tracked under dev/: %s\n' "$path" >&2
      status=1
      ;;
  esac

  bytes=$(wc -c < "$path")
  if (( bytes > max_source_bytes )); then
    printf 'Tracked dev/ file exceeds %d bytes: %s (%d bytes)\n' \
      "$max_source_bytes" "$path" "$bytes" >&2
    status=1
  fi

  if grep -Iq . "$path" && grep -Eq \
    '(/Users/[^/[:space:]]+|/home/[^/[:space:]]+|~/current_projects)' "$path"; then
    printf 'Personal absolute path is tracked under dev/: %s\n' "$path" >&2
    status=1
  fi
done < <(git ls-files -z -- dev/)

if (( status != 0 )); then
  printf 'dev/ source-only check failed. Put generated output under output/.\n' >&2
  exit "$status"
fi

printf 'dev/ source-only check passed.\n'
