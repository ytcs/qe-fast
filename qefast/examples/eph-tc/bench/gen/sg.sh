#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# DIAGNOSTIC: shared gauge. EPW only, on the off run's nscf and the off run's Wannier output (wannierize = .false.), with each set's
# ph.x data (save/). A control with its own save/ must give the off run's numbers. So any change left comes from the phonon/DFPT data, not from Wannier.
#  sg.sh <material> [sets...]
M=$1; shift; G=/root/gen; R=$G/run/$M; D=$G/diag/$M
source /root/qefast/env.sh; export QEF_PRESET=off
for S in off "$@"; do
  W=$D/sg_$S; rm -rf $W; mkdir -p $W; cd $W
  for f in $R/off/x.* $R/off/*.fmt $R/off/decay.*; do case "$f" in *.a2f*|*.dos|*.wout) ;; *) cp -r "$f" .;; esac; done
  [ -e $R/off/x.xml ] && cp $R/off/x.xml .
  cp -r $R/$S/save .
  sed "s/wannierize = .true./wannierize = .false./" $R/off/epw.in > epw.in
  rm -f restart.fmt selecq.fmt; mkdir -p x.ephmat
  nvidia-cuda-mps-control -d 2>/dev/null
  OMP_NUM_THREADS=1 mpirun --bind-to none -np 8 -x OMP_NUM_THREADS -x QEF_PRESET /root/qefast/qe-fast/build-gpu/bin/epw.x -npool 8 -in epw.in > epw.out 2>&1; rc=$?
  echo quit | nvidia-cuda-mps-control 2>/dev/null
  python3 - "$W" "$S" "$rc" <<'PY'
import sys
sys.path.insert(0, '/root/gen')
import qe_f2 as qe
W, n, rc = sys.argv[1:4]
try:
    m = qe.a2f_metrics(f"{W}/x.a2f_proj")
    print(f"sg {n:8s} rc {rc} lambda {m['lam']:.6f} wlog {m['wlog']:.5f} Tc {m['tc_ad']:.5f}")
except Exception as e:
    print("sg", n, "rc", rc, "no result", e)
PY
  rm -rf $W/x.save $W/save $W/*.epb* $W/x.epmatwp
done
echo SG_DONE $M
