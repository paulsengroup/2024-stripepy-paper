#!/usr/bin/env bash

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set -e
set -u
set -o pipefail

if ! git rev-parse --show-toplevel &> /dev/null; then
  1>&2 echo "This script must be run from within the repository root!"
  exit 1
fi

root="$(git rev-parse --show-toplevel)"
outdir="$root/benchmarks"
if compgen -G "$outdir/*.tsv" &> /dev/null; then
  1>&2 echo "Found one or more TSV files under folder \"$outdir\": please remove them before running this script."
  exit 1
fi

dataset_id='ENCFF993FGR'
resolution=10000
source_matrix="$root/data/$dataset_id.mcool"
dest_matrix="/tmp/$dataset_id.mcool"

if [ ! -f "$source_matrix" ]; then
  1>&2 echo "Unable to find matrix file \"$source_matrix\""
  exit 1
fi

mkdir -p "$outdir"

if [ "$(uname)" == 'Darwin' ]; then
  docker_user="$USER"
else
  docker_user='root'
fi

images=(
  stripepy:1.0.0
  chromosight:1.6.3
  stripenn:1.1.65.22
  stripecaller:0.1.0
)

for img in "${images[@]}"; do
  sudo -u "$docker_user" docker pull "ghcr.io/paulsengroup/2024-stripepy-paper/$img"
done

script="$(mktemp "${TMPDIR:-/tmp}/runner.XXXXXX.sh")"

# shellcheck disable=SC2064
trap "rm -f '$script'" EXIT

printf 'set -e\nset -u\nset -o pipefail\nset -x\n' >> "$script"

for img in "${images[@]}"; do
  tool_name="$(echo "$img" | cut -f 1 -d ':')"

  echo \
    docker run \
      --user "'$UID'" \
      -v "'$root/scripts/run_perf_benchmarks.py:/tmp/runme.py:ro'" \
      -v "'$source_matrix:$dest_matrix:ro'" \
      --entrypoint="'/opt/$tool_name/bin/python'" \
      --rm \
      "ghcr.io/paulsengroup/2024-stripepy-paper/$img" \
      /tmp/runme.py \
      "'$dest_matrix'" \
      "'$resolution'" \
      --tool "'$tool_name'" \|\
    tee "'$outdir/$tool_name.tsv'"

  echo "chown $UID '$outdir/$tool_name.tsv'"
done >> "$script"

cat "$script"

sudo -u "$docker_user" bash "$script"

scripts/summarize_perf_benchmarks.py "$outdir/"*.tsv
