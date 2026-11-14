#!/usr/bin/env python3
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
"""gamma_ref.py -- CPU reference for the impact-ionization rate of a hot electron (book E-13, issue #1).

  Im Sigma_nk = 1/(Nk Omega) sum_{q != 0} sum_{n' empty} phi^H Im W(q, w) phi,   w = e_nk - e_n'k-q > 0
  phi(G)      = <n, k| e^{i(q+G)r} |n', k-q>          (the orientation of chi0_ref's rho)
  Im W        = v^1/2 [(eps^-1 - eps^-1^H) / 2i] v^1/2  (RPA with local fields, chi0_ref.eps_inv)
  Gamma_ii    = 2 z |Im Sigma| / hbar,  here with z = 1

Every conduction state up to --emax-carrier above the conduction minimum is an initial state; the final states
are the empty states between the conduction minimum and the initial energy (holes are created in the valence
bands through Im W). Im W is built for each q on an energy grid and interpolated linearly. The q = 0 term
(head and wings of W) is left out here: the run reports the weight it leaves out, 1/Nk.

  gamma_ref.py <outdir> <prefix> [--emax-carrier eV] [--ecut-eps Ry] [--de eV] [--eta eV] [--scissor eV]
               [--out file] [--partial file.npz]
Writes n, k, E - E_cbm (eV), Im Sigma (eV), Gamma (1/s); prints Gamma averaged in 0.5 eV bins. --partial writes the
rate into every final state (n', k'), whose sum over final states is Gamma: the cascade draws the scattered carrier
from it.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chi0_ref import HA, System, binned_A, eps_inv   # noqa: E402

HBAR_HA = 2.4188843265857e-17   # s, hbar / E_h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir'); ap.add_argument('prefix')
    ap.add_argument('--emax-carrier', type=float, default=6.0, help='initial states up to this above E_cbm (eV)')
    ap.add_argument('--ecut-eps', type=float, default=4.0)
    ap.add_argument('--de', type=float, default=0.05)
    ap.add_argument('--eta', type=float, default=0.05)
    ap.add_argument('--out', default='gamma_ref.dat')
    ap.add_argument('--max-mem', type=float, default=2.0)
    ap.add_argument('--scissor', type=float, default=0.0, help='shift of the empty bands (eV), in W and the states')
    ap.add_argument('--partial', default='', help='also write the partial rates per final state (n\', k\') to this .npz')
    a = ap.parse_args()
    sy = System(a.outdir, a.prefix, scissor_ev=a.scissor)
    nv, nbnd, nk = sy.nv, sy.nbnd, sy.nk
    ecbm = sy.eig[:, nv].min() * HA
    evbm = sy.eig[:, nv - 1].max() * HA
    e = sy.eig * HA
    # initial states: conduction bands up to the carrier cutoff; the highest band used must lie above it everywhere
    nmax = nv
    while nmax < nbnd and e[:, nmax].min() <= ecbm + a.emax_carrier:
        nmax += 1
    if nmax >= nbnd - 1 or e[:, nbnd - 1].min() < ecbm + a.emax_carrier:
        sys.exit('gamma_ref: not enough bands above the carrier cutoff; raise nbnd in the nscf')
    wmax = a.emax_carrier + 0.5
    wgrid = np.arange(0.0, wmax + a.de, a.de)
    print(f'gamma_ref: {nk} k, bands {nv + 1}-{nmax} as initial states, E_cbm {ecbm:.3f} eV, gap {ecbm - evbm:.3f} eV, '
          f'W on {len(wgrid)} energies to {wmax} eV; q = 0 left out (weight {1 / nk:.4f})')
    imsig = np.zeros((nk, nbnd))
    part = []                                  # (n, k, n', k', Im Sigma contribution) for the cascade's final-state draw
    qlist = [iq for iq in range(nk) if np.abs(sy.kcr[iq] - np.round(sy.kcr[iq])).max() > 1e-8]
    for count, iq in enumerate(qlist):
        qcr = sy.kcr[iq] - np.ceil(sy.kcr[iq] - 0.5 - 1e-9)
        Gm, q2 = sy.gbasis(qcr, a.ecut_eps)
        NG = len(Gm)
        need = 16 * NG * NG * (int(np.ceil(wmax / a.de)) + 2 + 3 * len(wgrid)) / 1e9
        if need > a.max_mem:
            sys.exit(f'gamma_ref: {need:.2f} GB exceeds --max-mem {a.max_mem}')
        A = binned_A(sy, qcr, Gm, wmax, a.de, log=False)
        _, einv = eps_inv(sy, A, q2, a.de, wgrid, a.eta)
        sv = np.sqrt(4 * np.pi / q2)
        imW = sv[None, :, None] * ((einv - np.conj(np.swapaxes(einv, 1, 2))) / 2j) * sv[None, None, :]
        del einv, A
        kq, G0 = sy.kq_map(qcr)              # k' + q = k_kq + G0: final state at k', initial at k = k' + q
        for kp in range(nk):
            k = kq[kp]
            phi = sy.pairs(kp, k, G0[kp], slice(nv, nbnd), slice(nv, nmax), Gm)   # (n initial, n' final, G)
            for n in range(nv, nmax):
                w = e[k, n] - e[kp, nv:nbnd]                 # energy transfer to the pair
                ok = (w > 0) & (e[k, n] <= ecbm + a.emax_carrier)
                for jn in np.nonzero(ok)[0]:
                    x = w[jn] / a.de; j0 = int(x); f = x - j0
                    if j0 + 1 >= len(wgrid):
                        continue
                    M = (1 - f) * imW[j0] + f * imW[j0 + 1]
                    p = phi[n - nv, jn]
                    contrib = np.real(np.conj(p) @ M @ p)
                    imsig[k, n] += contrib
                    if a.partial:
                        part.append((n, k, nv + jn, kp, contrib))
        print(f'  q {count + 1}/{len(qlist)}  NG {NG}', flush=True)
    imsig /= nk * sy.omega
    gam = 2 * np.abs(imsig) / HBAR_HA
    if a.partial:
        pt = np.array(part)
        np.savez(a.partial, n=pt[:, 0].astype(int) + 1, k=pt[:, 1].astype(int) + 1, nf=pt[:, 2].astype(int) + 1,
                 kf=pt[:, 3].astype(int) + 1, gamma=2 * np.abs(pt[:, 4] / (nk * sy.omega)) / HBAR_HA,
                 e_init=e[pt[:, 1].astype(int), pt[:, 0].astype(int)] - ecbm,
                 e_final=e[pt[:, 3].astype(int), pt[:, 2].astype(int)] - ecbm)
    rows = [(n + 1, k + 1, e[k, n] - ecbm, imsig[k, n] * HA, gam[k, n])
            for k in range(nk) for n in range(nv, nmax) if e[k, n] <= ecbm + a.emax_carrier]
    np.savetxt(a.out, np.array(rows), fmt=['%d', '%d', '%.4f', '%.6e', '%.6e'],
               header='band k E-E_cbm(eV) ImSigma(eV) Gamma_ii(1/s), z = 1, q = 0 term left out')
    r = np.array(rows)
    print(f'gamma_ref: largest Im Sigma {r[:, 3].max():.2e} eV (must be <= 0)')
    for lo in np.arange(0.0, a.emax_carrier, 0.5):
        m = (r[:, 2] >= lo) & (r[:, 2] < lo + 0.5)
        if m.any():
            print(f'gamma_ref: E-E_cbm {lo:4.1f}-{lo + 0.5:4.1f} eV: {m.sum():4d} states, Gamma mean {r[m, 4].mean():.2e} 1/s, '
                  f'max {r[m, 4].max():.2e}')


if __name__ == '__main__':
    main()
