#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# CI gate for chi0.x: Si on a 4^3 grid without symmetry, 12 bands, Q = 0.5 (2pi/a, x); chi0.x against the
# reference qefast/gw/chi0_ref.py, every output column and the whole eps^-1 matrix, to 1e-8.
#   ci.sh <bin dir> <workdir>
set -euo pipefail
HERE=$(dirname "$(realpath "$0")")
BIN=$(realpath "$1"); mkdir -p "$2"; W=$(realpath "$2")
PSEUDO_DIR=$(realpath "$HERE/../../../../EPW/examples/sic/pp")
cd "$W"
sed -e "s|@PSEUDO@|$PSEUDO_DIR|" -e "s/8 8 8 0 0 0/4 4 4 0 0 0/" "$HERE/../si-eels/scf.in" > scf.in
sed -e "s/calculation = 'scf'/calculation = 'nscf'/" \
    -e "s/ecutwfc = 20.0/ecutwfc = 20.0\n    nbnd = 12\n    nosym = .true.\n    noinv = .true./" \
    -e "s/conv_thr = 1.0d-10/conv_thr = 1.0d-10\n    diago_full_acc = .true./" scf.in > nscf.in
"$BIN/pw.x" -in scf.in > scf.out
"$BIN/pw.x" -in nscf.in > nscf.out
cat > chi0.in <<EOI
&inputchi0
  prefix = 'si', outdir = './out', q = 0.5, 0.0, 0.0, ecut_eps = 4.0, emax = 20.0, de = 0.05, eta = 0.1
  fileout = 'chi0-f90.dat', filfull = 'einv-f90.bin'
/
EOI
"$BIN/chi0.x" -in chi0.in > chi0.out
grep 'chi0:' chi0.out
python3 "$HERE/../../../gw/chi0_ref.py" out si --q 0.5 0 0 --ecut-eps 4 --emax 20 --de 0.05 --eta 0.1 \
  --out chi0-py.dat --full einv-py.npy | grep 'chi0_ref:'
python3 - <<'EOP'
import numpy as np, sys
a, b = np.loadtxt('chi0-py.dat'), np.loadtxt('chi0-f90.dat')
d_col = np.abs(a - b).max()
e = np.load('einv-py.npy'); nw, ng, _ = e.shape
f = np.swapaxes(np.fromfile('einv-f90.bin', dtype=np.complex128).reshape(nw, ng, ng), 1, 2)
rd = lambda p: [l.split()[1:4] for l in open(p) if l.strip() and l[0] != '#']
same_g = rd('einv-py.npy.G.txt') == rd('einv-f90.bin.G.txt')
d_full = np.abs(e - f).max()
print(f'chi0 gate: columns max |diff| {d_col:.2e}, eps^-1 ({nw} x {ng}^2) max |diff| {d_full:.2e}, same G list {same_g}')
sys.exit(0 if (d_col < 1e-8 and d_full < 1e-8 and same_g) else 1)
EOP
