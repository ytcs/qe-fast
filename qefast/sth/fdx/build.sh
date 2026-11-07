#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# Build the cuFFTDx microbenchmark: qefast/sth/fdx/build.sh [sm, e.g. 86 or 90] -> build/fdxbench
# Needs NVHPC 26.5 (nvfortran, nvcc), MathDx (pip nvidia-mathdx 25.6) and CUTLASS headers:
#   MDX=<.../nvidia/mathdx/include> CUTLASS=<cutlass/include>
set -euo pipefail
cd "$(dirname "$(realpath "$0")")"
SM=${1:-86}
NV=/opt/nvidia/hpc_sdk/Linux_x86_64/26.5
export PATH=$NV/cuda/13.2/bin:$NV/compilers/bin:$PATH
MDX=${MDX:-/root/mdx/lib/python3.12/site-packages/nvidia/mathdx/include}
CUTLASS=${CUTLASS:-/root/cutlass/include}
B=build${EMU:+-emu}; mkdir -p $B && cd $B
nvcc -std=c++17 -O3 -arch=sm_$SM -DFDX_SM=${SM}0 ${EMU:+-DFDX_EMU} -I$MDX -I$CUTLASS -Xptxas -v -c ../sth_fftdx.cu -o sth_fftdx.o 2> ptxas.txt || { cat ptxas.txt; exit 1; }
grep -E "spill" ptxas.txt | grep -v " 0 bytes spill stores, 0 bytes spill loads" || true
nvcc -std=c++17 -O3 -arch=sm_$SM ${EMU:+-DFDX_EMU} -I$NV/math_libs/13.2/include -c ../fdxbench.cu -o fdxbench.o
nvcc -arch=sm_$SM fdxbench.o sth_fftdx.o -o fdxbench -L$NV/math_libs/13.2/lib64 -Xlinker -rpath=$NV/math_libs/13.2/lib64 -lcufft
echo "built $(pwd)/fdxbench"
nvcc -O3 -arch=sm_$SM -I$NV/math_libs/13.2/include ../gemmbench.cu -o gemmbench -L$NV/math_libs/13.2/lib64 -Xlinker -rpath=$NV/math_libs/13.2/lib64 -lcublas && echo "built $(pwd)/gemmbench"
