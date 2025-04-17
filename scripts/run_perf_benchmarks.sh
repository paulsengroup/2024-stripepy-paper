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
if ! compgen -G "$outdir/*.tsv"; then
  1>&2 echo "Found one or more TSV files under folder $outdir: please remove them before running this script."
  exit 1
fi

mkdir -p "$outdir"

if [ "$(uname)" == "Darwin" ]; then
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


sudo -u "$docker_user" run \
  --user "$UID" \
  -v "$root/scripts/run_perf_benchmarks.py:/tmp/runme.py:ro" \
  -v "$root/data/4DNFI9GMP2J8.mcool:/tmp/4DNFI9GMP2J8.mcool:ro" \
  --entrypoint=/opt/stripepy/bin/python \
  --rm \
  ghcr.io/paulsengroup/2024-stripepy-paper/stripepy:1.0.0 \
  /tmp/runme.py \
  /tmp/4DNFI9GMP2J8.mcool \
  10000 \
  --suppress-output \
  --tool stripepy |
  tee "$outdir/stripepy.tsv"


sudo -u "$docker_user" run \
  --user "$UID" \
  -v "$root/scripts/run_perf_benchmarks.py:/tmp/runme.py:ro" \
  -v "$root/data/4DNFI9GMP2J8.mcool:/tmp/4DNFI9GMP2J8.mcool:ro" \
  --entrypoint=/opt/chromosight/bin/python \
  --rm \
  ghcr.io/paulsengroup/2024-stripepy-paper/chromosight:1.6.3 \
  /tmp/runme.py \
  /tmp/4DNFI9GMP2J8.mcool \
  10000 \
  --suppress-output \
  --tool chromosight |
  tee "$outdir/chromosight.tsv"

sudo -u "$docker_user" run \
  --user "$UID" \
  -v "$root/scripts/run_perf_benchmarks.py:/tmp/runme.py:ro" \
  -v "$root/data/4DNFI9GMP2J8.mcool:/tmp/4DNFI9GMP2J8.mcool:ro" \
  --entrypoint=/opt/stripenn/bin/python \
  --rm \
  ghcr.io/paulsengroup/2024-stripepy-paper/stripenn:1.1.65.22 \
  /tmp/runme.py \
  /tmp/4DNFI9GMP2J8.mcool \
  10000 \
  --suppress-output \
  --tool stripenn |
  tee "$outdir/stripenn.tsv"

sudo -u "$docker_user" run \
  --user "$UID" \
  -v "$root/scripts/run_perf_benchmarks.py:/tmp/runme.py:ro" \
  -v "$root/data/4DNFI9GMP2J8.mcool:/tmp/4DNFI9GMP2J8.mcool:ro" \
  --entrypoint=/opt/stripecaller/bin/python \
  --rm \
  ghcr.io/paulsengroup/2024-stripepy-paper/stripecaller:0.1.0 \
  /tmp/runme.py \
  /tmp/4DNFI9GMP2J8.mcool \
  10000 \
  --suppress-output \
  --tool stripecaller |
  tee "$outdir/stripecaller.tsv"
