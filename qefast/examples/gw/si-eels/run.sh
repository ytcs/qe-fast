#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# Si loss function -Im 1/eps(Q, w) with turbo_eels (Liouville-Lanczos), one run per Q.
#   [APPROX=TDDFT|IPA|RPA_with_CLFE] run.sh <bin dir> <workdir> [Q in 2pi/a along x ...]   (default Q: 0.1 2.1)
# Prints, per Q, the loss peak and turbo_spectrum's f-sum violation; spectra go to <workdir>/Q<q>/si.plot_eps.dat.
set -euo pipefail
HERE=$(dirname "$(realpath "$0")")
BIN=$(realpath "$1"); mkdir -p "$2"; W=$(realpath "$2"); shift 2
QS=("${@:-0.1 2.1}"); [ $# -eq 0 ] && QS=(0.1 2.1)
PSEUDO_DIR=$(realpath "$HERE/../../../../EPW/examples/sic/pp")     # Si.pz-vbc.UPF, in the tree
cd "$W"
sed "s|@PSEUDO@|$PSEUDO_DIR|" "$HERE/scf.in" > scf.in
"$BIN/pw.x" -in scf.in > scf.out
for q in "${QS[@]}"; do
  d=Q$q; rm -rf "$d"; mkdir "$d"; cp -r out "$d/"
  ( cd "$d"
    sed "s/@ITERMAX@/1500/;s/@Q1@/$q/;s/@Q2@/0.0/;s/@Q3@/0.0/;s/@APPROX@/${APPROX:-TDDFT}/" "$HERE/eels.in" > eels.in
    sed "s/@ITERMAX@/1500/;s/@EMAX@/150.0/" "$HERE/spectrum.in" > spectrum.in
    "$BIN/turbo_eels.x" -in eels.in > eels.out
    "$BIN/turbo_spectrum.x" -in spectrum.in > spectrum.out
    peak=$(awk '!/#/ && $3 > m {m = $3; e = $1} END {printf "%.2f eV (%.3f)", e, m}' si.plot_eps.dat)
    echo "Q = $q 2pi/a: loss peak $peak; $(grep -o 'violation.*' spectrum.out)" )
done
