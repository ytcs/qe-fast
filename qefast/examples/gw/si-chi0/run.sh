#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# Si inverse dielectric matrix by sum over bands (qefast/gw/chi0_ref.py): scf, nscf on the full 8^3 grid without
# symmetry, then chi0_ref.py per Q. Same Si as ../si-eels, so the turbo_eels runs there are the reference.
#   [NBND=60] [ECUT_EPS=8] run.sh <bin dir> <workdir> [Q in 2pi/a along x ...]   (default Q: 0.25 2.25)
# Q must be commensurate with the 8^3 grid (0.25 and 2.25 are; 0.1 and 2.1 are not).
set -euo pipefail
HERE=$(dirname "$(realpath "$0")")
BIN=$(realpath "$1"); mkdir -p "$2"; W=$(realpath "$2"); shift 2
QS=("$@"); [ $# -eq 0 ] && QS=(0.25 2.25)
PSEUDO_DIR=$(realpath "$HERE/../../../../EPW/examples/sic/pp")
cd "$W"
sed "s|@PSEUDO@|$PSEUDO_DIR|" "$HERE/../si-eels/scf.in" > scf.in
"$BIN/pw.x" -in scf.in > scf.out
sed -e "s/calculation = 'scf'/calculation = 'nscf'/" \
    -e "s/ecutwfc = 20.0/ecutwfc = 20.0\n    nbnd = ${NBND:-60}\n    nosym = .true.\n    noinv = .true./" \
    -e "s/conv_thr = 1.0d-10/conv_thr = 1.0d-10\n    diago_full_acc = .true./" scf.in > nscf.in
"$BIN/pw.x" -in nscf.in > nscf.out
for q in "${QS[@]}"; do
  python3 "$HERE/../../../gw/chi0_ref.py" out si --q "$q" 0 0 --ecut-eps "${ECUT_EPS:-8}" --emax 80 --de 0.02 \
    --eta 0.1 --out "chi0-Q$q.dat" | grep -v '^  '
done
