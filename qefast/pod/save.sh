#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# Save the dev tree (source, build-gpu objects and binaries, pseudo, NVTX shim) to the build cache.
# Run at the end of a dev session so the next pod starts warm. Runs and data are not saved.
set -euo pipefail
CACHE=${CACHE:-/cache}; CC=${CC:-80}; ARC=$CACHE/qefast-sm$CC.tar.zst
cd /root/qefast
tar -I 'zstd -3 -T0' -cf "$ARC.tmp" --exclude='qe-fast/build-cpu' qe-fast nvtxshim fdxdeps env.sh
mv "$ARC.tmp" "$ARC" && echo "saved $(du -h "$ARC" | cut -f1) to $ARC"
