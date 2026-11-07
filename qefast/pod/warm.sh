#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# Warm start for a dev pod.
# Expects the build-cache volume mounted at $CACHE (default /cache).
#   curl -sSfL https://raw.githubusercontent.com/ytcs/qe-fast/main/qefast/pod/warm.sh | bash
# Cache hit: restore source + build-gpu (NVTX on) + pseudo + shim, git pull, incremental rebuild.
# Cache miss: full setup + build, then save to the cache.
set -euo pipefail
CACHE=${CACHE:-/cache}; J=${J:-28}; CC=${CC:-80}; ARC=$CACHE/qefast-sm$CC.tar.zst
NV=/opt/nvidia/hpc_sdk/Linux_x86_64/26.5
apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq zstd curl git libfftw3-dev pkg-config cmake ninja-build python3 python3-numpy >/dev/null
mkdir -p /root/qefast && cd /root/qefast
cat > env.sh <<EOI
NV=$NV
export PATH=\$NV/comm_libs/mpi/bin:\$NV/compilers/bin:\$NV/profilers/13.2/Nsight_Systems/bin:\$PATH
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export LIBRARY_PATH=/root/qefast/nvtxshim:\${LIBRARY_PATH:-} LD_LIBRARY_PATH=/root/qefast/nvtxshim:\${LD_LIBRARY_PATH:-}
EOI
if [ -f "$ARC" ]; then
  echo "warm: restoring $(du -h "$ARC" | cut -f1) from cache"
  tar -I 'zstd -d -T0' -xf "$ARC" -C /root/qefast
  cd qe-fast && git fetch -q origin && git checkout -q "${BRANCH:-main}" && git pull -q --ff-only origin "${BRANCH:-main}" && cd ..
else
  echo "warm: cache miss, full setup"
  [ -d qe-fast ] || git clone -q -b "${BRANCH:-main}" https://github.com/ytcs/qe-fast.git
  (cd qe-fast && git submodule update -q --init external/wannier90 external/mbd external/fox)
  qe-fast/qefast/pseudo.sh
fi
source env.sh
eval "$(FDX_DEPS=/root/qefast/fdxdeps qe-fast/qefast/sth/fdx/deps.sh)"   # without these the FFTDX path is the stub
NVTX=1 CPU=0 J=$J MDX=$MDX CUTLASS=$CUTLASS qe-fast/qefast/pod/build-gpu.sh qe-fast $CC
[ -f "$ARC" ] || CC=$CC qe-fast/qefast/pod/save.sh
echo WARM_OK
