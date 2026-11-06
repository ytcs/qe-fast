#!/usr/bin/env python3
"""wfcentres.py <run dir> [materials...]: Cartesian Wannier centres of each preset's run against the material's `off` run.

The spread sum cannot see a WF move to another lattice cell, and lambda depends on those home cells, so identity of two
Wannier solutions is checked on the Cartesian centres (no lattice reduction). Centres are matched one to one
(Hungarian); the report gives the largest matched distance in Angstrom, how many WFs sit further apart than TOL, and the
spread sums. Same-index distance is also given, since a run with the same gauge keeps the order.
"""
import glob
import os
import re
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

TOL = 0.02   # Angstrom
PAT = re.compile(r"WF centre and spread\s+\d+\s+\(\s*(\S+?),\s*(\S+?),\s*(\S+?)\s*\)\s+(\S+)")


def final_state(d):
    """Centres and spreads of the last 'Final State' block of x.wout."""
    t = open(f"{d}/x.wout").read()
    t = t[t.rindex("Final State"):]
    t = t[:t.index("Sum of centres and spreads")]
    v = np.array([[float(x) for x in m] for m in PAT.findall(t)])
    return v[:, :3], v[:, 3]


def main():
    r = sys.argv[1]
    mats = sys.argv[2:] or sorted(os.path.basename(d) for d in glob.glob(f"{r}/G*"))
    for m in mats:
        if not os.path.exists(f"{r}/{m}/off/x.wout"):
            print(f"{m}: no off x.wout")
            continue
        c0, s0 = final_state(f"{r}/{m}/off")
        for p in ("safe", "fast", "max"):
            if not os.path.exists(f"{r}/{m}/{p}/x.wout"):
                continue
            c1, s1 = final_state(f"{r}/{m}/{p}")
            dist = np.linalg.norm(c0[:, None, :] - c1[None, :, :], axis=2)
            i, j = linear_sum_assignment(dist)
            d = dist[i, j]
            same = np.linalg.norm(c0 - c1, axis=1)
            print(f"{m} {p}: nWF {len(d)} matched max {d.max():.6f} A, over {TOL} A: {(d > TOL).sum()}; same-index max {same.max():.6f} A; "
                  f"spread sum {s0.sum():.6f} -> {s1.sum():.6f}, per-WF spread max diff (matched) {np.abs(s0[i] - s1[j]).max():.6f}")


if __name__ == "__main__":
    main()
