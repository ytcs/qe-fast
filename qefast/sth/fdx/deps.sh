#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# MathDx (cuFFTDx) and CUTLASS headers for the fftdx engine path. Prints the two include dirs.
#   eval "$(qefast/sth/fdx/deps.sh)"   -> MDX=..., CUTLASS=...
set -uo pipefail
D=${FDX_DEPS:-/root/qefast/fdxdeps}
mkdir -p $D
if ! ls $D/mdx/lib/python3*/site-packages/nvidia/mathdx/include/cufftdx.hpp > /dev/null 2>&1; then
  python3 -m venv $D/mdx >&2 && $D/mdx/bin/pip install -q nvidia-mathdx==25.6.0 >&2 \
    || echo "deps.sh: MathDx install failed (python3-venv missing?); the build falls back to the FFTDX stub" >&2
fi
[ -d $D/cutlass/include ] || git clone -q --depth 1 -b v3.9.2 https://github.com/NVIDIA/cutlass.git $D/cutlass >&2
M=$(ls -d $D/mdx/lib/python3*/site-packages/nvidia/mathdx/include 2>/dev/null || true)
if [ -z "$M" ]; then
  echo "deps.sh: no MathDx headers: fused FFT kernels (FFTDX) will be off, B1's ph stage is then ~2.8x slower" >&2
  exit 1
fi
echo "export MDX=$M CUTLASS=$D/cutlass/include"
