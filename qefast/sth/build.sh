#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# Build the offline harness of the batched Sternheimer engine (no CMake, no QE libraries).
#   qefast/sth/build.sh            -> qefast/sth/build/sth_harness.x
# Needs nvfortran (NVHPC) with CUDA; on the dev pod: source /root/qefast/env.sh first.
set -euo pipefail
cd "$(dirname "$(realpath "$0")")"
FC=${FC:-nvfortran}
FLAGS=${FLAGS:-"-O2 -cuda -acc -gpu=cc80 -cudalib=cufft,cublas"}
mkdir -p build
$FC $FLAGS -module build -c ../config/qef_config.f90 -o build/qef_config.o && $FC $FLAGS -module build -c sth_engine.f90 -o build/sth_engine.o
$FC $FLAGS -module build -c sth_harness.f90 -o build/sth_harness.o
nvc -O2 -c fdx/sth_fftdx_stub.c -o build/sth_fftdx_stub.o
$FC $FLAGS build/qef_config.o build/sth_engine.o build/sth_harness.o build/sth_fftdx_stub.o -o build/sth_harness.x
echo "built $(pwd)/build/sth_harness.x"
