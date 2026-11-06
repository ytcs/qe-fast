#!/usr/bin/env python3
"""Held-out generalisation check: scf/nscf/ph inputs of G1-G8 and D1-D3.
  mkinputs.py <out dir> <archive dir (F0c cells)> [B13C2 relaxed scf.in]
epw.in needs E_F and the bands: mkepw.py writes it after the preparation nscf."""
import json, math, os, re, sys
import numpy as np
import qe_f2 as qe

OUT, ARC = sys.argv[1], sys.argv[2]
B13 = sys.argv[3] if len(sys.argv) > 3 else None


def fcc(a):
    return [[-a / 2, 0, a / 2], [0, a / 2, a / 2], [-a / 2, a / 2, 0]]


def bcc(a):
    return [[a / 2, a / 2, a / 2], [-a / 2, a / 2, a / 2], [-a / 2, -a / 2, a / 2]]


def cell(path):
    t = open(path).read()
    lat = [[float(x) for x in l.split()] for l in
           re.search(r"CELL_PARAMETERS angstrom\n((?:.*\n){3})", t).group(1).strip().split("\n")]
    m = re.search(r"ATOMIC_POSITIONS crystal\n(.*?)K_POINTS", t, re.S).group(1).strip().split("\n")
    return dict(lattice=lat, species=[l.split()[0] for l in m], frac=[[float(x) for x in l.split()[1:4]] for l in m])


def b13c2(a=5.19, alpha=65.7, R=1.70, u=0.385):
    """(B12)(CBC), R-3m: a regular icosahedron on the rhombohedral vertex with six polar atoms along +-a_i (the 5-fold axes
    near the cell edges), six equatorial atoms, and the C-B-C chain along the 3-fold axis. A starting geometry: vc-relaxed."""
    ca = math.cos(math.radians(alpha))
    cb = math.sqrt((1 + 2 * ca) / 3)            # cos of the angle between a_i and the 3-fold axis
    sb = math.sqrt(1 - cb * cb)
    lat = [[a * sb * math.cos(p), a * sb * math.sin(p), a * cb] for p in (0, 2 * math.pi / 3, 4 * math.pi / 3)]
    L = np.array(lat)
    t1, t2 = math.radians(37.3774), math.radians(79.1877)
    pts = []
    for i in range(3):
        p = 2 * math.pi * i / 3
        pts.append((math.sin(t1) * math.cos(p), math.sin(t1) * math.sin(p), math.cos(t1)))
        p2 = p + math.pi / 3
        pts.append((math.sin(t2) * math.cos(p2), math.sin(t2) * math.sin(p2), math.cos(t2)))
    pts = [R * np.array(v) for v in pts]
    pts += [-v for v in pts]
    frac = [list(np.linalg.solve(L.T, v) % 1.0) for v in pts]
    sp = ["B"] * 12 + ["B", "C", "C"]
    frac += [[0.5, 0.5, 0.5], [u, u, u], [1 - u, 1 - u, 1 - u]]
    return dict(lattice=lat, species=sp, frac=[[float(x) for x in f] for f in frac])


M = {
    "G1": dict(st=dict(lattice=fcc(4.25), species=["Mg", "O"], frac=[[0, 0, 0], [.5, .5, .5]]), k=(6, 6, 6), q=(2, 2, 2), nk=(4, 4, 4), epw=None),
    "G2": dict(st=dict(lattice=fcc(4.04), species=["Al"], frac=[[0, 0, 0]]), k=(12, 12, 12), q=(4, 4, 4), nk=(8, 8, 8), kf=(16,) * 3, qf=(8,) * 3),
    "G3": dict(st=dict(lattice=bcc(3.31), species=["Nb"], frac=[[0, 0, 0]]), k=(12, 12, 12), q=(4, 4, 4), nk=(8, 8, 8), kf=(16,) * 3, qf=(8,) * 3),
    "G4": dict(st=cell(f"{ARC}/H3S/scf.in"), k=(15, 15, 15), q=(3, 3, 3), nk=(6, 6, 6), kf=(16,) * 3, qf=(8,) * 3, ecut=72),
    "G5": dict(st=cell(f"{ARC}/LiAl2Rh/scf.in"), k=(8, 8, 8), q=(2, 2, 2), nk=(4, 4, 4), kf=(10,) * 3, qf=(5,) * 3, ecut=88),
    "G6": dict(st=dict(lattice=fcc(4.25), species=["Ti", "N"], frac=[[0, 0, 0], [.5, .5, .5]]), k=(8, 8, 8), q=(2, 2, 2), nk=(4, 4, 4), kf=(12,) * 3, qf=(6,) * 3,
               sys="nr1 = 39, nr2 = 39, nr3 = 39"),
    "G7": dict(st=cell(f"{ARC}/MgB3C3__s1/scf.in"), k=(4, 4, 3), q=(2, 2, 1), nk=(4, 4, 2), kf=(8,) * 3, qf=(4,) * 3, ecut=84),
    "G8": dict(st=cell(B13) if B13 else b13c2(), k=(4, 4, 4), q=(2, 2, 2), nk=(4, 4, 4), kf=(6,) * 3, qf=(3,) * 3),
    "D1": dict(st=dict(lattice=bcc(2.87), species=["Fe"], frac=[[0, 0, 0]]), k=(6, 6, 6), q=(2, 2, 2), nspin=2, ph_only=True),
    "D2": dict(st=dict(lattice=fcc(4.25), species=["Mg", "O"], frac=[[0, 0, 0], [.5, .5, .5]]), k=(6, 6, 6), q=(2, 2, 2), fixed=True, ph_only=True),
    "D3": dict(st=cell(f"{ARC}/MgB3C3__s1/scf.in"), k=(6, 6, 7), q=(2, 2, 2), ecut=84, ph_only=True),
}


def kblock(n):
    pts = [(i / n[0], j / n[1], l / n[2]) for i in range(n[0]) for j in range(n[1]) for l in range(n[2])]
    return f"K_POINTS crystal\n{len(pts)}\n" + "".join(f" {a:.10f} {b:.10f} {c:.10f} {1 / len(pts):.6e}\n" for a, b, c in pts)


def fix(t, m):
    if m.get("fixed"):
        t = re.sub(r"occupations = 'smearing'\n\s*smearing = 'mp'\n\s*degauss = [\d.]+\n", "occupations = 'fixed'\n", t)
    return t


for name, m in M.items():
    d = f"{OUT}/{name}"
    os.makedirs(d, exist_ok=True)
    st = m["st"]
    ecut = m.get("ecut") or qe.ecut_rule(st["species"])
    es = m.get("sys", "")
    open(f"{d}/scf.in", "w").write(fix(qe.pw_input(st, "scf", "x", m["k"], ecut, extra_sys=es, nspin=m.get("nspin", 1)), m))
    masses = "".join(f"  amass({i + 1}) = {qe.MASS[e]}\n" for i, e in enumerate(qe.species_order(st)))
    q = m["q"]
    last = "  start_q = 1\n  last_q = 2\n" if m.get("ph_only") else ""
    open(f"{d}/ph.in", "w").write(f"&inputph\n  prefix = 'x'\n  fildyn = 'x.dyn'\n{masses}  outdir = './'\n  ldisp = .true.\n"
                                  f"  fildvscf = 'dvscf'\n  nq1 = {q[0]}, nq2 = {q[1]}, nq3 = {q[2]}\n  tr2_ph = 1.0d-14\n  reduce_io = .true.\n{last}/\n")
    if m.get("ph_only"):
        continue
    nw = qe.n_wann(st["species"])
    nbnd = None
    if m.get("epw", True) is not None:
        ns = qe.n_shell(st["species"], os.environ["PSEUDO_DIR"])
        nbnd = ns + nw + max(8, nw // 2)
    else:
        nbnd = 16
    open(f"{d}/nscf.in", "w").write(qe.pw_input(st, "nscf", "x", None, ecut, nbnd=nbnd, extra_sys=es, kblock=kblock(m["nk"])))
    json.dump(dict(species=st["species"], q=m["q"], nk=m["nk"], kf=m.get("kf"), qf=m.get("qf"), epw=m.get("epw", True) is not None, nw=nw), open(f"{d}/meta.json", "w"))
    print(name, "atoms", len(st["species"]), "ecut", ecut, "nw", nw, "nbnd", nbnd)
