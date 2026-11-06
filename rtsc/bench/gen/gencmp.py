#!/usr/bin/env python3
"""gencmp.py <run dir> [materials...]: the pre-registered gate and report of groups/compute/GEN-LOG.md. For every material and set against
that material's all-off run: frequencies (Gamma acoustic excluded), lambda, omega_log, Allen-Dynes Tc (a2f_proj, funnel a2f_metrics), margin used,
stage speedups, peak memory, decline lines, Wannier spread sum."""
import glob, json, os, re, statistics as st, sys
import numpy as np
import qe_f2 as qe

R = sys.argv[1]
MATS = sys.argv[2:] or sorted(os.path.basename(d) for d in glob.glob(f"{R}/G*"))
SETS = ["safe", "fast", "max"]


def freqs(d):
    out = {}
    for f in glob.glob(f"{d}/x.dyn*"):
        m = re.search(r"\.dyn(\d+)$", f)
        if not m or m.group(1) == "0":
            continue
        t = open(f).read()
        q = re.search(r"q = \(\s*(\S+)\s+(\S+)\s+(\S+)", t)
        w = np.array([float(x) for x in re.findall(r"freq.*=\s*\S+\s*\[THz\]\s*=\s*(\S+)\s*\[cm-1\]", t)])
        out[int(m.group(1))] = (np.array([float(x) for x in q.groups()]), w)
    return out


def times(d):
    t = {}
    p = f"{d}/timings.tsv"
    if os.path.exists(p):
        for l in open(p):
            a = l.split("\t")
            if len(a) == 2:
                t[a[0]] = float(a[1])
    return t


def a2f(d):
    p = f"{d}/x.a2f_proj"
    if not os.path.exists(p):
        p = next(iter(glob.glob(f"{d}/x.a2f*")), None)
    if not p:
        return None
    try:
        return qe.a2f_metrics(p)
    except Exception as e:
        return dict(err=str(e))


def spread(d):
    p = f"{d}/x.wout"
    if os.path.exists(p):
        v = re.findall(r"Sum of centres and spreads.*?([\d.]+)\s*$", open(p).read(), re.M)
        if v:
            return float(v[-1])
    return None


def lines(d):
    o = []
    for f in ("ph.out", "nscf.out", "scf.out", "epw.out"):
        p = f"{d}/{f}"
        if not os.path.exists(p):
            continue
        seen = set()
        for l in open(p, errors="replace"):
            if re.search(r"rtsc cbsr auto|unsupported grid|engine declined|CG kept|sth_engine: (needs|irrep needs)|Error in routine|rtsc nscf symrot: rep|vanilla", l) \
                    and "q = Gamma" not in l:
                k = re.sub(r"[\d.]+", "#", l.strip())[:60]
                if k not in seen:
                    seen.add(k)
                    o.append(f"{f}: {l.strip()[:230]}")
    return o


rows = []
sp = {}
for M in MATS:
    A = f"{R}/{M}/alloff"
    ra = open(f"{A}/status").read().strip() if os.path.exists(f"{A}/status") else "?"
    fa, ta, aa, sa = freqs(A), times(A), a2f(A), spread(A)
    print(f"\n## {M}: all-off rc {ra}, total {ta.get('total', 0):.1f} s, peak {open(A + '/peakmem').read().strip() if os.path.exists(A + '/peakmem') else '?'} MiB, nq {len(fa)}, "
          f"lambda {aa and aa.get('lam')}, wlog {aa and aa.get('wlog')}, Tc {aa and aa.get('tc_ad')}, spread {sa}")
    print(" ", "stages s:", {k: round(v, 1) for k, v in ta.items()})
    for l in lines(A):
        print("   all-off:", l)
    for S in SETS:
        D = f"{R}/{M}/{S}"
        if not os.path.exists(f"{D}/status"):
            print(f"  {S}: not run")
            continue
        rc = open(f"{D}/status").read().strip()
        ft, tt, at, stt = freqs(D), times(D), a2f(D), spread(D)
        worst, wn, ok = 0.0, "", rc == "0"
        dmax = 0.0
        if set(fa) != set(ft):
            ok = False
            wn = f"q sets differ {len(fa)} vs {len(ft)}"
        for i in sorted(set(fa) & set(ft)):
            xq, a = fa[i]
            b = ft[i][1]
            mask = np.ones(len(a), bool)
            if np.allclose(xq, 0):
                mask[np.argsort(np.abs(a))[:3]] = False
            tol = np.maximum(0.005 * np.abs(a), 1.0)
            r = np.abs(a - b) / tol
            r[~mask] = 0
            j = int(np.argmax(r))
            dmax = max(dmax, np.abs(a - b)[mask].max() if mask.any() else 0)
            if r[j] > worst:
                worst = float(r[j])
                wn = f"q{i} mode {j + 1} ({a[j]:.1f} cm-1, |dw| {abs(a[j] - b[j]):.3f})"
        res = {}
        if aa and at and "lam" in aa and "lam" in at:
            for k, name, rel, ab in (("lam", "lambda", .01, 0), ("wlog", "wlog", .01, 0), ("tc_ad", "Tc", .02, .5)):
                tol = max(rel * abs(aa[k]), ab) or 1e-12
                x = abs(at[k] - aa[k]) / tol
                res[name] = (at[k], (at[k] - aa[k]) / aa[k] if aa[k] else 0.0)
                if x > worst:
                    worst, wn = x, f"{name} {aa[k]:.4f} -> {at[k]:.4f}"
        elif M != "G1":
            ok = False
            wn += " (no a2f)"
        ok = ok and worst <= 1.0
        spd = {k: (ta[k] / tt[k] if k in tt and tt[k] > 0 and k in ta else float("nan")) for k in ("scf", "ph", "nscf", "epw", "total")}
        pk = open(f"{D}/peakmem").read().strip() if os.path.exists(f"{D}/peakmem") else "?"
        print(f"  {S}: rc {rc} gate {'PASS' if ok else 'FAIL'} margin {100 * worst:.1f} % [{wn}] max|dw| {dmax:.3f}; speedup " +
              " ".join(f"{k} {v:.2f}" for k, v in spd.items()) + f"; peak {pk} MiB; spread {stt}; " +
              " ".join(f"{k} {v[1]:+.2e}" for k, v in res.items()))
        for l in lines(D):
            print("     ", l)
        sp.setdefault(S, []).append((spd["total"], M, ok, spd["ph"]))
print("\n## summary")
for S in SETS:
    v = [x for x in sp.get(S, []) if x[0] == x[0]]
    if v:
        tot = [x[0] for x in v]
        print(f"{S}: n {len(v)} median {st.median(tot):.2f} worst {min(tot):.2f} ({min(v)[1]}) best {max(tot):.2f}; gate pass {sum(x[2] for x in v)}/{len(v)}; "
              f"ph median {st.median([x[3] for x in v]):.2f}")
