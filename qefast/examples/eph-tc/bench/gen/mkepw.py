#!/usr/bin/env python3
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
"""mkepw.py <input dir> <prep dir>: epw.in from the preparation scf + nscf (funnel setup A: atomic projections, G1, frozen window)."""
import json, os, re, sys
import qe_f2 as qe
d, p = sys.argv[1], sys.argv[2]
m = json.load(open(f"{d}/meta.json"))
xs = qe.read_xml(f"{p}/scf/x.xml")
xn = qe.read_xml(f"{p}/nscf/x.xml")
xn["ef"] = xs["ef"]
sp = m["species"]
ns, info = qe.n_semicore(xs, sp, os.environ["PSEUDO_DIR"])
nw = m["nw"]
w = qe.frozen_window(xn, ns, nw)
if w is None:
    raise SystemExit(f"{d}: no frozen window with <= {nw} bands")
t = qe.epw_input("x", sp, 1, 1, 1, 1, ns, nw, w, xs["ef"])
n, q, kf, qf = m["nk"], m["q"], m["kf"], m["qf"]
t = t.replace("nk1 = 1, nk2 = 1, nk3 = 1", f"nk1 = {n[0]}, nk2 = {n[1]}, nk3 = {n[2]}")
t = t.replace("nq1 = 1, nq2 = 1, nq3 = 1", f"nq1 = {q[0]}, nq2 = {q[1]}, nq3 = {q[2]}")
t = t.replace("nkf1 = 1, nkf2 = 1, nkf3 = 1", f"nkf1 = {kf[0]}, nkf2 = {kf[1]}, nkf3 = {kf[2]}")
t = t.replace("nqf1 = 1, nqf2 = 1, nqf3 = 1", f"nqf1 = {qf[0]}, nqf2 = {qf[1]}, nqf3 = {qf[2]}")
open(f"{d}/epw.in", "w").write(t)
print(os.path.basename(d), "E_F", round(xs["ef"], 3), "nsemi", ns, "nw", nw, "window +", w, info)
