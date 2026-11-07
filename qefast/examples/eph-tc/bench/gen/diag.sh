#!/usr/bin/env bash
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
# DIAGNOSTIC: which lever moves the Wannier solution. Chains on one material, all-off plus one lever group.
#  diag.sh <material>
M=$1; G=/root/gen; Q=/root/qefast/qe-fast; D=$G/diag/$M
source /root/qefast/env.sh
one() { name=$1; shift
  W=$D/$name; rm -rf $W; mkdir -p $W
  ( export QEF_PRESET=off; for e in "$@"; do export "$e"; done; [ "$M" = G1 ] && export QEF_RUN_NOEPW=1
    flock /root/qefast/TIMING.lock timeout -k 20 2700 bash $Q/qefast/examples/eph-tc/bench/run-chain.sh $G/in/$M $Q/build-gpu/bin $W 8 8 > $W/run.log 2>&1 ); rc=$?
  python3 - "$W" "$name" "$rc" <<'PY'
import re, sys
sys.path.insert(0, '/root/gen')
import qe_f2 as qe
W, n, rc = sys.argv[1:4]
try:
    v = re.findall(r"Sum of centres and spreads.*?([\d.]+)\s*$", open(f"{W}/x.wout").read(), re.M)[-1]
    m = qe.a2f_metrics(f"{W}/x.a2f_proj")
    print(f"{n:14s} rc {rc} spread {float(v):.6f} lambda {m['lam']:.5f} wlog {m['wlog']:.4f} Tc {m['tc_ad']:.4f}")
except Exception as e:
    print(n, "rc", rc, "no result", e)
PY
  rm -rf $W/_ph0 $W/x.epmatwp $W/*.epb* $W/*.ephmat*
}
one alloff2
one alloff3
one pwsymrot QEF_PW_SYMROT=1
one w90 QEF_W90_NT=8 QEF_W90_WNT=8 QEF_W90_MMN=1 QEF_W90_AMN=1 QEF_W90_NOHR=1
one w90nt QEF_W90_NT=8
one w90wnt QEF_W90_WNT=8
one w90mmn QEF_W90_MMN=1
one w90amn QEF_W90_AMN=1
one epw QEF_EPW_FDX=1 QEF_EPW_BEC=1 QEF_EPW_STAR=1 QEF_EPW_B2W=1 QEF_EPW_IDLEWAIT=1 QEF_EPW_NOVME=1
one ph QEF_STH_FFTDX=1 QEF_STH_FDXDP=1 QEF_STH_OBATCH=1 QEF_STH_QCACHE=1 QEF_NSCF_SYMROT=1 QEF_PHQ_FUSE=1 QEF_PHQ_VKB=1 QEF_STH_DVB=1 QEF_DVXC_GPU=1 QEF_STH_KDEV=2 QEF_STH_PSIR_PART=1 QEF_STH_DRHO2=1
echo DIAG_DONE $M
