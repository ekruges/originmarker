#!/usr/bin/env bash
# Run the OriginMarker CLI on tesla, against the array corpus that lives there.
#
# WHY TESLA. 12 cores and 20 GB free against the NAS's 4 and roughly 1, and node 22 is already
# installed at /opt/node22. The NAS runs the live site and its memory headroom is the binding
# constraint there; here it is not a constraint at all.
#
# WHY THE JUMP HOST. tesla is wifi-only behind a NAT guest bridge, so this Mac cannot route to
# 192.168.1.40 directly even when tesla is up. newton sits on both networks and reaches it in
# 1.5 ms, so ~/.ssh/config sends tesla through it. Without that, a live machine looks like a dead
# one: the failure is a connection timeout either way.
#
# THE CORPUS IS ALREADY THERE, at /opt/om/data, uncompressed, from an earlier transfer in this
# project. 884 arrays: DIETER 93, JENNA 178, ROBLES 264, TREFF 349. Check before copying anything;
# a second copy was most of the way to being made before anyone looked.
#
#   usage: scripts/om-tesla.sh <om subcommand and flags...>
#     scripts/om-tesla.sh census /opt/om/data
#     scripts/om-tesla.sh cohort /opt/om/data/ROBLES --ref /root/parent.probes --json
set -euo pipefail
REPO=/root/originmarker
HEAP="${OM_HEAP:-12288}"
ssh tesla "cd $REPO && /opt/node22/bin/node --experimental-strip-types \
  --max-old-space-size=$HEAP cli/om.ts $(printf '%q ' "$@")"
