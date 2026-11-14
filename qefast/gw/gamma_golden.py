# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check of the impact-ionisation contraction (gamma_ref.contract, and so gamma_ii.x) against the golden rule.

To lowest order in the bare Coulomb interaction, Im W = v Im chi0 v. Then the self-energy route
    Gamma = 2 |Im Sigma| / hbar,   Im Sigma = 1/(Nk Omega) sum_{q, n'} contract(phi, Im W)
must equal Fermi's golden rule for the decay into an electron (n', k') and an electron-hole pair (c k1+q, v k1):
    Gamma = (2 pi / hbar) * 2 * sum |M|^2 delta(w - D),   M = 1/V sum_G v(q+G) <n'k'|e^{-i(q+G)r}|nk> <c k1+q|e^{i(q+G)r}|v k1>
(V = Nk Omega, 2 for the spin of the pair, delta the same binned triangle as chi0, interpolated in w the same way).
This fixes the prefactors and the order of the G indices: Im W_GG' is not symmetric, and phi^H Im W phi instead of
phi^T Im W phi^* was 26 % off on Si 4^3.

  gamma_golden.py <outdir> <prefix> [--emax-carrier 2.5] [--ecut-eps 4] [--de 0.05] [--eta 0.05]; exit 1 on failure
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chi0_ref import HA, System, binned_A, kk_kernel   # noqa: E402
from gamma_ref import HBAR_HA, WMIN, contract          # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('outdir'); ap.add_argument('prefix')
ap.add_argument('--emax-carrier', type=float, default=2.5)
ap.add_argument('--ecut-eps', type=float, default=4.0)
ap.add_argument('--de', type=float, default=0.05)
ap.add_argument('--eta', type=float, default=0.05)
a = ap.parse_args()
sy = System(a.outdir, a.prefix)
nv, nbnd, nk = sy.nv, sy.nbnd, sy.nk
e = sy.eig * HA
ecbm = e[:, nv].min()
nmax = nv
while nmax < nbnd and e[:, nmax].min() <= ecbm + a.emax_carrier:
    nmax += 1
wmax = a.emax_carrier + 0.5
wgrid = np.arange(0.0, wmax + a.de, a.de)
nb = int(np.ceil(wmax / a.de)) + 2
pref = sy.spin / (nk * sy.omega)
# Im[K(w_i - D_b) - K(w_i + D_b)] for every bin b on the grid (1/Hartree), the kernel of chi0
Kb = np.array([np.imag(kk_kernel((wgrid - b * a.de) / HA, a.eta / HA, True)
                       - kk_kernel((wgrid + b * a.de) / HA, a.eta / HA, True)) for b in range(nb)])
se = np.zeros((nk, nbnd))
gr = np.zeros((nk, nbnd))
for iq in range(nk):
    if np.abs(sy.kcr[iq] - np.round(sy.kcr[iq])).max() < 1e-8:
        continue
    qcr = sy.kcr[iq] - np.ceil(sy.kcr[iq] - 0.5 - 1e-9)
    Gm, q2 = sy.gbasis(qcr, a.ecut_eps)
    v = 4 * np.pi / q2
    A = binned_A(sy, qcr, Gm, wmax, a.de, log=False)
    imW0 = v[None, :, None] * (pref * np.einsum('bi,bgh->igh', Kb, A)) * v[None, None, :]
    kq, G0 = sy.kq_map(qcr)
    rho, dw = [], []                                   # the transitions of chi0 and their delta weights on the grid
    for ik in range(nk):
        r = sy.pairs(ik, kq[ik], G0[ik], slice(0, nv), slice(nv, nbnd), Gm)
        D = (sy.eig[kq[ik], nv:nbnd][:, None] - sy.eig[ik, :nv][None, :]) * HA / a.de
        for ic in range(nbnd - nv):
            for iv in range(nv):
                j0 = int(np.floor(D[ic, iv])); f = D[ic, iv] - j0
                if j0 + 1 < nb:
                    rho.append(r[ic, iv]); dw.append(-((1 - f) * Kb[j0] + f * Kb[j0 + 1]) / np.pi)
    rho, dw = np.array(rho), np.array(dw)
    for kp in range(nk):
        k = kq[kp]
        phi = sy.pairs(kp, k, G0[kp], slice(nv, nbnd), slice(nv, nmax), Gm)
        for n in range(nv, nmax):
            if e[k, n] > ecbm + a.emax_carrier:
                continue
            for jn in range(nbnd - nv):
                w = e[k, n] - e[kp, nv + jn]
                x = w / a.de; j0 = int(x); f = x - j0
                if w <= WMIN or j0 + 1 >= len(wgrid):
                    continue
                p = phi[n - nv, jn]
                se[k, n] += contract(p, (1 - f) * imW0[j0] + f * imW0[j0 + 1])
                if len(rho):
                    gr[k, n] += np.sum(np.abs(rho @ (v * np.conj(p))) ** 2 * ((1 - f) * dw[:, j0] + f * dw[:, j0 + 1]))
V = nk * sy.omega
g_se = 2 * np.abs(se / V) / HBAR_HA
g_gr = 2 * np.pi / HBAR_HA * 2 * gr / V ** 2
m = g_gr > 1e-9 * g_gr.max()
d = np.abs(g_se[m] / g_gr[m] - 1).max() if m.any() else np.inf
print(f'gamma_golden: {m.sum()} states, self-energy route / golden rule - 1: max {d:.2e} '
      f'(rates {g_gr[m].min():.2e} .. {g_gr[m].max():.2e} 1/s)')
sys.exit(0 if (m.any() and d < 1e-8) else 1)
