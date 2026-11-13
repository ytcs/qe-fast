#!/usr/bin/env python3
# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
"""chi0_ref.py -- CPU reference for the full inverse dielectric matrix eps^-1_GG'(q, w) (issue #1).

Sum over bands on a pw.x nscf run (norm-conserving, full k grid without symmetry, k + q on the grid):

  rho_t(G)      = <n', k+q| e^{i(q+G)r} |n, k>                      pair densities, any band pair
  A_GG'(D)      = sum_t rho_t(G)^* rho_t(G') delta(D - D_t)           binned in D (Hermitian, PSD per bin)
  chi0_GG'(w)   = 2/(Nk Omega) sum_bins A [1/(w - D + i eta) - 1/(w + D + i eta)]   (time-reversal symmetric)
  eps_GG'       = delta_GG' - v(q+G)^1/2 chi0_GG' v(q+G')^1/2,  v = 4 pi / |q+G|^2
  loss L(Q, w)  = -Im eps^-1_{GQ GQ},  Q = q + GQ;  without local fields -Im 1/(1 - v chi0_{GQ GQ}).
  turbo_eels' 'IPA' drops the Hartree term, so its loss is -Im(1 + v chi0) = Im eps without local fields.

Hartree atomic units inside; energies in and out in eV. This is the reference that the GPU program is checked
against; it is written for clarity, not speed.

  chi0_ref.py <outdir> <prefix> --q Qx Qy Qz [--ecut-eps Ry] [--emax eV] [--de eV] [--eta eV] [--out file]
Q in units of 2 pi / alat, Cartesian. Writes w, L_IPA, L_RPA, Re eps_M (RPA) and prints the loss peak and the
f-sum fraction int w L dw / (pi/2 w_p^2) with w_p^2 = 4 pi n_valence.
"""
import argparse
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

HA = 27.211386245988  # eV


def read_xml(save):
    r = ET.parse(os.path.join(save, 'data-file-schema.xml')).getroot()
    out = r.find('output')
    alat = float(out.find('atomic_structure').get('alat'))
    cell = np.array([[float(x) for x in out.find('atomic_structure/cell/' + a).text.split()] for a in ('a1', 'a2', 'a3')])
    nr = [int(out.find('basis_set/fft_grid').get(k)) for k in ('nr1', 'nr2', 'nr3')]
    bs = out.find('band_structure')
    nelec = float(bs.find('nelec').text)
    nbnd = int(bs.find('nbnd').text) if bs.find('nbnd') is not None else int(bs.find('nbnd_up').text)
    noncolin = bs.find('noncolin').text.strip() == 'true'
    ks = bs.findall('ks_energies')
    xk = np.array([[float(x) for x in k.find('k_point').text.split()] for k in ks])   # Cartesian, 2 pi / alat
    eig = np.array([[float(x) for x in k.find('eigenvalues').text.split()] for k in ks])  # Hartree
    # decline anything but norm-conserving pseudopotentials
    for sp in out.findall('atomic_species/species'):
        pf = sp.find('pseudo_file').text.strip()
        txt = open(os.path.join(save, pf), errors='replace').read(20000) if os.path.exists(os.path.join(save, pf)) else ''
        if 'pseudo_type="US"' in txt or 'pseudo_type="PAW"' in txt or 'PP_PAW' in txt:
            sys.exit(f'chi0_ref: {pf} is ultrasoft or PAW; only norm-conserving pseudopotentials are supported')
    return dict(alat=alat, cell=cell, nr=nr, nelec=nelec, nbnd=nbnd, noncolin=noncolin, xk=xk, eig=eig)


def read_wfc(path):
    """QE binary wfc file: (ik, xk, ispin, gamma_only, scalef), (ngw, igwx, npol, nbnd), b1-b3, mill, evc(:, n)."""
    with open(path, 'rb') as f:
        def rec():
            n = np.frombuffer(f.read(4), np.int32)[0]
            b = f.read(n)
            f.read(4)
            return b
        h = rec()
        gamma_only = np.frombuffer(h[32:36], np.int32)[0] != 0
        if gamma_only:
            sys.exit('chi0_ref: gamma_only wavefunctions are not supported')
        ngw, igwx, npol, nbnd = np.frombuffer(rec(), np.int32)
        rec()
        mill = np.frombuffer(rec(), np.int32).reshape(igwx, 3)
        evc = np.array([np.frombuffer(rec(), np.complex128) for _ in range(nbnd)])   # (nbnd, npol * igwx)
    return mill, evc.reshape(nbnd, npol, igwx), npol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir'); ap.add_argument('prefix')
    ap.add_argument('--q', nargs=3, type=float, required=True)
    ap.add_argument('--ecut-eps', type=float, default=4.0, help='|q+G|^2 cutoff of the local-field basis (Ry)')
    ap.add_argument('--emax', type=float, default=60.0); ap.add_argument('--de', type=float, default=0.02)
    ap.add_argument('--eta', type=float, default=0.1); ap.add_argument('--out', default='chi0_ref.dat')
    ap.add_argument('--nbnd', type=int, default=0, help='use only the lowest nbnd bands (convergence tests)')
    ap.add_argument('--max-mem', type=float, default=2.0, help='refuse to run above this estimate (GB)')
    a = ap.parse_args()
    save = os.path.join(a.outdir, a.prefix + '.save')
    d = read_xml(save)
    alat, cell, nr = d['alat'], d['cell'], np.array(d['nr'])
    omega = abs(np.linalg.det(cell))
    bg = 2 * np.pi * np.linalg.inv(cell).T          # rows b_i, bohr^-1
    nk, nbnd = d['eig'].shape
    if a.nbnd:
        nbnd = min(nbnd, a.nbnd)
    nv = int(round(d['nelec'] / (1 if d['noncolin'] else 2)))
    kcr = d['xk'] * (2 * np.pi / alat) @ cell.T / (2 * np.pi)   # crystal coordinates of k
    Qcart = np.array(a.q) * 2 * np.pi / alat
    Qcr = Qcart @ cell.T / (2 * np.pi)
    GQ = np.ceil(Qcr - 0.5 - 1e-9).astype(int)      # q = Q - GQ in (-1/2, 1/2] (crystal)
    qcr = Qcr - GQ
    # k + q on the grid
    kq = np.empty(nk, int); G0 = np.empty((nk, 3), int)
    for ik in range(nk):
        t = kcr[ik] + qcr
        dd = t[None, :] - kcr
        j = np.argmin(np.abs(dd - np.round(dd)).sum(1))
        if np.abs(dd[j] - np.round(dd[j])).max() > 1e-6:
            sys.exit('chi0_ref: k + q is not on the k grid; choose q commensurate with the grid')
        kq[ik] = j; G0[ik] = np.round(dd[j]).astype(int)
    # local-field basis: G with |q+G|^2 <= ecut_eps (Ry), q+G in bohr^-1; GQ is in it by construction if ecut allows
    n = int(np.ceil(np.sqrt(a.ecut_eps) / np.min(np.linalg.norm(bg, axis=1)))) + 2
    rng = np.arange(-n, n + 1)
    Gm = np.array(np.meshgrid(rng, rng, rng, indexing='ij')).reshape(3, -1).T
    qG = (qcr + Gm) @ bg
    q2 = (qG ** 2).sum(1)
    keep = q2 <= a.ecut_eps            # |q+G|^2 in bohr^-2 is the kinetic energy in Ry
    keep |= np.all(Gm == GQ, axis=1)
    Gm, q2 = Gm[keep], q2[keep]
    order = np.argsort(q2); Gm, q2 = Gm[order], q2[order]
    iQ = int(np.where(np.all(Gm == GQ, axis=1))[0][0])
    NG = len(Gm)
    v = 4 * np.pi / q2
    print(f'chi0_ref: {nk} k, {nbnd} bands ({nv} occupied), FFT {nr.tolist()}, {NG} G (|q+G|^2 <= {a.ecut_eps} Ry), '
          f'Q crystal {np.round(Qcr, 4).tolist()} = q {np.round(qcr, 4).tolist()} + G {GQ.tolist()}')
    # bins; memory: the binned matrices plus one omega chunk of chi / eps / eps^-1, estimated before allocating
    nb = int(np.ceil(a.emax / a.de)) + 2
    nw = len(np.arange(0.0, a.emax, a.de))
    wc = max(1, min(nw, int(2.5e8 // (16 * NG * NG * 4 + 16 * nb))))       # ~0.25 GB per chunk
    need = (16 * nb * NG * NG + wc * (16 * NG * NG * 4 + 16 * nb)) / 1e9
    print(f'chi0_ref: memory estimate {need:.2f} GB ({nb} bins x {NG}^2, omega chunks of {wc})')
    if need > a.max_mem:
        sys.exit(f'chi0_ref: {need:.2f} GB exceeds --max-mem {a.max_mem} GB; raise --de, lower --emax or --ecut-eps')
    A = np.zeros((nb, NG, NG), complex)
    gidx = tuple((Gm % nr).T)
    for ik in range(nk):
        mk, ck, npol = read_wfc(os.path.join(save, f'wfc{ik + 1}.dat'))
        mq, cq, _ = read_wfc(os.path.join(save, f'wfc{kq[ik] + 1}.dat'))
        # u_{n,k}(r) and u_{n',k+q}(r); the k+q coefficients sit at mill - G0
        def to_r(mill, c, nb_):
            box = np.zeros((nb_, npol) + tuple(nr), complex)
            ix = tuple((mill % nr).T)
            for s in range(npol):
                box[(slice(None), s) + ix] = c[:nb_, s, :]
            return np.fft.ifftn(box, axes=(2, 3, 4)) * np.prod(nr)
        uk = to_r(mk, ck[:nv], nv)
        uq = to_r(mq - G0[ik], cq[nv:nbnd], nbnd - nv)
        ek, eq = d['eig'][ik], d['eig'][kq[ik]]
        for iv in range(nv):
            f = (np.conj(uq) * uk[iv][None]).sum(1)                   # (nc, nr1, nr2, nr3)
            rho = np.fft.ifftn(f, axes=(1, 2, 3))[(slice(None),) + gidx]   # (nc, NG) = <c,k+q|e^{i(q+G)r}|v,k>
            D = (eq[nv:nbnd] - ek[iv]) * HA
            x = D / a.de; j0 = np.floor(x).astype(int); w1 = x - j0
            ok = j0 + 1 < nb
            for jj, ww in ((j0, 1 - w1), (j0 + 1, w1)):
                for b in np.unique(jj[ok]):
                    m = ok & (jj == b)
                    A[b] += (np.conj(rho[m]) * ww[m][:, None]).T @ rho[m]   # sum_t rho_t(G)^* rho_t(G')
        if (ik + 1) % 64 == 0:
            print(f'  {ik + 1}/{nk} k', flush=True)
    # Hermiticity / positivity of each bin (free check)
    herm = max(np.abs(A[b] - A[b].conj().T).max() for b in range(nb))
    mineig = min(np.linalg.eigvalsh(0.5 * (A[b] + A[b].conj().T)).min() for b in range(nb))
    print(f'chi0_ref: max |A - A^H| {herm:.2e}, min eigenvalue of the bins {mineig:.2e}')
    Db = np.arange(nb) * a.de / HA
    w = np.arange(0.0, a.emax, a.de) / HA
    eta = a.eta / HA
    pref = 2.0 / (nk * omega)
    Af = A.reshape(nb, -1)
    sv = np.sqrt(v)
    eQQ = np.empty(len(w), complex); einvQQ = np.empty(len(w), complex)
    for i0 in range(0, len(w), wc):
        ww = w[i0:i0 + wc]
        K = 1.0 / (ww[:, None] - Db[None, :] + 1j * eta) - 1.0 / (ww[:, None] + Db[None, :] + 1j * eta)
        chi = pref * (K @ Af).reshape(len(ww), NG, NG)   # time reversal: the antiresonant term has the same A
        eps = np.eye(NG)[None] - sv[None, :, None] * chi * sv[None, None, :]
        eQQ[i0:i0 + wc] = eps[:, iQ, iQ]
        einvQQ[i0:i0 + wc] = np.linalg.inv(eps)[:, iQ, iQ]
    L_rpa = -einvQQ.imag
    L_ipa = -(1.0 / eQQ).imag
    im_eps0 = eQQ.imag      # = -v Im chi0_QQ: what turbo_eels calls the IPA loss (no Hartree term)
    epsM = 1.0 / einvQQ
    wev = w * HA
    wp2 = 4 * np.pi * d['nelec'] / omega
    np.savetxt(a.out, np.column_stack([wev, im_eps0, L_ipa, L_rpa, epsM.real, epsM.imag]),
               header='w(eV)  Im eps_NLF(=turbo IPA)  L_NLF  L_RPA  Re eps_M(RPA)  Im eps_M(RPA)')
    print(f'chi0_ref: Im eps without local fields (turbo IPA): peak {wev[np.argmax(im_eps0)]:.2f} eV ({im_eps0.max():.3f}), '
          f'f-sum fraction {np.trapezoid(w * im_eps0, w) / (0.5 * np.pi * wp2):.3f}')
    for name, L in (('NLF', L_ipa), ('RPA', L_rpa)):
        fs = np.trapezoid(w * L, w) / (0.5 * np.pi * wp2)
        print(f'chi0_ref: {name} loss peak {wev[np.argmax(L)]:.2f} eV ({L.max():.3f}), f-sum fraction {fs:.3f}')
    print(f'chi0_ref: Re eps_M(w=0) RPA {epsM[0].real:.3f}, no local fields {eQQ[0].real:.3f}')


if __name__ == '__main__':
    main()
