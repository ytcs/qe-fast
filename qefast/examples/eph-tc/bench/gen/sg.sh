#!/usr/bin/env bash
# DIAGNOSTIC (labelled; GEN-LOG verdict rule): shared gauge. EPW only, on all-off's nscf and all-off's Wannier output (wannierize = .false.), with each set's
# ph.x data (save/). A control with all-off's own save/ must give all-off's numbers. So any change left comes from the phonon/DFPT data, not from Wannier.
#  sg.sh <material> [sets...]
M=$1; shift; G=/root/gen; R=$G/run/$M; D=$G/diag/$M
source /root/qefast/env.sh; source $G/env/alloff.env
for S in alloff "$@"; do
  W=$D/sg_$S; rm -rf $W; mkdir -p $W; cd $W
  for f in $R/alloff/x.* $R/alloff/*.fmt $R/alloff/decay.*; do case "$f" in *.a2f*|*.dos|*.wout) ;; *) cp -r "$f" .;; esac; done
  [ -e $R/alloff/x.xml ] && cp $R/alloff/x.xml .
  cp -r $R/$S/save .
  sed "s/wannierize = .true./wannierize = .false./" $R/alloff/epw.in > epw.in
  nvidia-cuda-mps-control -d 2>/dev/null
  OMP_NUM_THREADS=1 mpirun --bind-to none -np 8 -x OMP_NUM_THREADS -x QEF_STH -x QEF_STRICT -x QEF_W90_NT /root/qefast/bin-gen/epw.x -npool 8 -in epw.in > epw.out 2>&1; rc=$?
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
