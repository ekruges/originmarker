#!/usr/bin/env bash
# Run the OriginMarker CLI on the NAS, against the array corpus that lives there.
#
# WHY NOT JUST RUN IT LOCALLY. The 884-array corpus is 36 GB and used to live on a USB stick, which
# means it is absent whenever the stick is not plugged in and invisible to any process that cannot
# see that volume. It now lives on the NAS at /volume1/docker/originmarker-data, which is always up.
#
# WHY A CONTAINER. The NAS host has node v18 and this CLI needs 22 for --experimental-strip-types.
# A node:22-slim container gives exactly the right runtime and changes nothing about the host,
# which matters on a box also serving the live site.
#
# The repo it runs is the one every deploy already copies to /volume1/docker/originmarker, so this
# is the same code the web tool serves rather than a second checkout that can drift.
#
#   usage: scripts/om-nas.sh <om subcommand and flags...>
#     scripts/om-nas.sh census /data
#     scripts/om-nas.sh stage "/data/SNP array data/DIETER/<file>.probes"
#     scripts/om-nas.sh cohort "/data/SNP array data/ROBLES" --ref /data/parent.probes --json
#
# Paths inside the container: /w is the repo, /data is the corpus. Memory is the one real limit,
# so --stride 8 for anything that sweeps the whole corpus; the NAS has about 1 GB free.
set -euo pipefail
REPO=/volume1/docker/originmarker
DATA=/volume1/docker/originmarker-data
MEM="${OM_NAS_HEAP:-3072}"

ssh nas "cd $REPO && docker run --rm \
  -v $REPO:/w -v $DATA:/data -w /w \
  node:22-slim node --experimental-strip-types --max-old-space-size=$MEM \
  cli/om.ts $(printf '%q ' "$@")"
