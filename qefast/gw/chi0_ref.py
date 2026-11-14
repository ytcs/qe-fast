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
trapz = getattr(np, 'trapezoid', None) or np.trapz   # numpy < 2 (Ubuntu 24.04's python3-numpy)


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
    # supported: fixed occupations (insulators), no spin polarisation, no magnetisation (time reversal is assumed)
    occ = bs.find('occupations_kind')
    if occ is None or occ.text.strip() != 'fixed':
        sys.exit('chi0_ref: only occupations = fixed (insulators) are supported')
    if bs.find('lsda').text.strip() == 'true':
        sys.exit('chi0_ref: spin-polarised (lsda) runs are not supported')
    mag = out.find('magnetization/absolute')
    if mag is not None and abs(float(mag.text)) > 1e-6:
        sys.exit('chi0_ref: magnetic runs are not supported (the antiresonant term assumes time reversal)')
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


class System:
    """A pw.x nscf run: cell, k list (crystal), eigenvalues (Hartree), FFT grid, occupied-band count."""

    def __init__(self, outdir, prefix, nbnd=0):
        self.save = os.path.join(outdir, prefix + '.save')
        d = read_xml(self.save)
        self.__dict__.update(d)
        self.nr = np.array(d['nr'])
        self.omega = abs(np.linalg.det(self.cell))
        self.bg = 2 * np.pi * np.linalg.inv(self.cell).T          # rows b_i, bohr^-1
        self.nk, nb = self.eig.shape
        self.nbnd = min(nb, nbnd) if nbnd else nb
        self.nv = int(round(self.nelec / (1 if self.noncolin else 2)))
        self.spin = 1.0 if self.noncolin else 2.0                # spin degeneracy only for collinear wavefunctions
        self.kcr = self.xk * (2 * np.pi / self.alat) @ self.cell.T / (2 * np.pi)

    def to_crystal(self, Qcart_2pi_alat):
        return np.array(Qcart_2pi_alat) * 2 * np.pi / self.alat @ self.cell.T / (2 * np.pi)

    def kq_map(self, qcr):
        """k + q = k_j + G0 for every k: (j, G0)."""
        kq = np.empty(self.nk, int); G0 = np.empty((self.nk, 3), int)
        for ik in range(self.nk):
            dd = (self.kcr[ik] + qcr)[None, :] - self.kcr
            j = np.argmin(np.abs(dd - np.round(dd)).sum(1))
            if np.abs(dd[j] - np.round(dd[j])).max() > 1e-6:
                sys.exit('chi0_ref: k + q is not on the k grid; choose q commensurate with the grid')
            kq[ik] = j; G0[ik] = np.round(dd[j]).astype(int)
        return kq, G0

    def gbasis(self, qcr, ecut, must=None):
        """G (crystal) with |q+G|^2 <= ecut (bohr^-2 = Ry), sorted by |q+G|; returns G, |q+G|^2."""
        n = int(np.ceil(np.sqrt(ecut) / np.min(np.linalg.norm(self.bg, axis=1)))) + 2
        rng = np.arange(-n, n + 1)
        Gm = np.array(np.meshgrid(rng, rng, rng, indexing='ij')).reshape(3, -1).T
        q2 = (((qcr + Gm) @ self.bg) ** 2).sum(1)
        keep = q2 <= ecut
        if must is not None:
            keep |= np.all(Gm == must, axis=1)
        Gm, q2 = Gm[keep], q2[keep]
        o = np.argsort(np.round(q2 * 1e10), kind='stable')   # ties broken by enumeration order, as chi0.x does
        return Gm[o], q2[o]

    def u_r(self, ik, bands, shift=None):
        """Periodic parts u_n(r) of bands (slice) at k-point ik, on the FFT grid; coefficients moved by -shift."""
        mill, c, npol = read_wfc(os.path.join(self.save, f'wfc{ik + 1}.dat'))
        if shift is not None:
            mill = mill - shift
        c = c[bands]
        box = np.zeros((c.shape[0], npol) + tuple(self.nr), complex)
        ix = tuple((mill % self.nr).T)
        for sp in range(npol):
            box[(slice(None), sp) + ix] = c[:, sp, :]
        return np.fft.ifftn(box, axes=(2, 3, 4)) * np.prod(self.nr)

    def pairs(self, ik, jk, G0, ket, bra, Gm):
        """rho[b, a, G] = <bra_b, k+q| e^{i(q+G)r} |ket_a, k>, k+q = k_jk + G0 (ket at ik, bra at jk)."""
        uk = self.u_r(ik, ket)
        uq = self.u_r(jk, bra, shift=G0)
        gidx = tuple((Gm % self.nr).T)
        out = np.empty((uq.shape[0], uk.shape[0], len(Gm)), complex)
        for ia in range(uk.shape[0]):
            f = (np.conj(uq) * uk[ia][None]).sum(1)
            out[:, ia, :] = np.fft.ifftn(f, axes=(1, 2, 3))[(slice(None),) + gidx]
        return out


def binned_A(sy, qcr, Gm, emax, de, log=True):
    """A[b] = sum over valence -> conduction transitions at transfer q of rho^* rho^T, linear binning in D (eV)."""
    nb = int(np.ceil(emax / de)) + 2
    NG = len(Gm)
    A = np.zeros((nb, NG, NG), complex)
    kq, G0 = sy.kq_map(qcr)
    nv, nbnd = sy.nv, sy.nbnd
    for ik in range(sy.nk):
        rho = sy.pairs(ik, kq[ik], G0[ik], slice(0, nv), slice(nv, nbnd), Gm)     # (nc, nv, NG)
        ek, eq = sy.eig[ik], sy.eig[kq[ik]]
        for iv in range(nv):
            r = rho[:, iv, :]
            x = (eq[nv:nbnd] - ek[iv]) * HA / de
            j0 = np.floor(x).astype(int); w1 = x - j0
            ok = j0 + 1 < nb
            for jj, ww in ((j0, 1 - w1), (j0 + 1, w1)):
                for b in np.unique(jj[ok]):
                    m = ok & (jj == b)
                    A[b] += (np.conj(r[m]) * ww[m][:, None]).T @ r[m]          # sum_t rho_t(G)^* rho_t(G')
        if log and (ik + 1) % 64 == 0:
            print(f'  {ik + 1}/{sy.nk} k', flush=True)
    return A


def eps_inv(sy, A, q2, de, w_ev, eta_ev):
    """eps and eps^-1 (local fields, RPA) at the energies w_ev (eV); returns (eps, eps^-1), each (nw, NG, NG)."""
    nb = A.shape[0]
    Db = np.arange(nb) * de / HA
    w = np.asarray(w_ev) / HA
    eta = eta_ev / HA
    K = 1.0 / (w[:, None] - Db[None, :] + 1j * eta) - 1.0 / (w[:, None] + Db[None, :] + 1j * eta)
    chi = (sy.spin / (sy.nk * sy.omega)) * (K @ A.reshape(nb, -1)).reshape(len(w), *A.shape[1:])  # same A both terms
    sv = np.sqrt(4 * np.pi / q2)
    eps = np.eye(len(q2))[None] - sv[None, :, None] * chi * sv[None, None, :]
    return eps, np.linalg.inv(eps)


def chi_kk_fft(sy, A, de, nw, eta_ev):
    """chi0 on w_i = i de (i < nw) by FFT: both terms depend on w - D or w + D only, so on the uniform grid they are a
    discrete convolution and a correlation of A_b with K(m) = 1/(m de + i eta). Returns chi (nw, NG, NG); equal to
    eps_inv's direct sum up to rounding, at O(N_G^2 N log N) instead of O(N_G^2 N_w N_b)."""
    nb = A.shape[0]
    Af = A.reshape(nb, -1)
    m1 = np.arange(-(nb - 1), nw)                       # w_i - D_b = (i - b) de
    m2 = np.arange(0, nw + nb - 1)                      # w_i + D_b = (i + b) de
    K1 = 1.0 / (m1 * de / HA + 1j * eta_ev / HA)
    K2 = 1.0 / (m2 * de / HA + 1j * eta_ev / HA)
    L = 1 << int(np.ceil(np.log2(2 * nb + nw)))
    F1 = np.fft.fft(K1, L)[:, None]
    F2 = np.fft.fft(K2, L)[:, None]
    chi = np.empty((nw, Af.shape[1]), complex)
    block = max(1, int(4e6 // L))                       # FFT temporaries <= ~64 MB
    for c0 in range(0, Af.shape[1], block):
        a = Af[:, c0:c0 + block]
        t1 = np.fft.ifft(np.fft.fft(a, L, axis=0) * F1, axis=0)[nb - 1:nb - 1 + nw]
        t2 = np.fft.ifft(np.fft.fft(a[::-1], L, axis=0) * F2, axis=0)[nb - 1:nb - 1 + nw]
        chi[:, c0:c0 + block] = t1 - t2
    return (sy.spin / (sy.nk * sy.omega)) * chi.reshape(nw, *A.shape[1:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir'); ap.add_argument('prefix')
    ap.add_argument('--q', nargs=3, type=float, required=True)
    ap.add_argument('--ecut-eps', type=float, default=4.0, help='|q+G|^2 cutoff of the local-field basis (Ry)')
    ap.add_argument('--emax', type=float, default=60.0); ap.add_argument('--de', type=float, default=0.02)
    ap.add_argument('--eta', type=float, default=0.1); ap.add_argument('--out', default='chi0_ref.dat')
    ap.add_argument('--nbnd', type=int, default=0, help='use only the lowest nbnd bands (convergence tests)')
    ap.add_argument('--max-mem', type=float, default=2.0, help='refuse to run above this estimate (GB)')
    ap.add_argument('--full', default='', help='also write the whole eps^-1_GG\'(w) to this .npy (memory-mapped) '
                    'and the G list to <name>.G.txt')
    a = ap.parse_args()
    sy = System(a.outdir, a.prefix, a.nbnd)
    Qcr = sy.to_crystal(a.q)
    GQ = np.ceil(Qcr - 0.5 - 1e-9).astype(int)      # q = Q - GQ in (-1/2, 1/2] (crystal)
    qcr = Qcr - GQ
    if np.abs(qcr).max() < 1e-8:
        sys.exit('chi0_ref: q = 0 (Q on the reciprocal lattice) needs the head and wings treatment; not supported')
    Gm, q2 = sy.gbasis(qcr, a.ecut_eps, must=GQ)
    iQ = int(np.where(np.all(Gm == GQ, axis=1))[0][0])
    NG = len(Gm)
    print(f'chi0_ref: {sy.nk} k, {sy.nbnd} bands ({sy.nv} occupied), FFT {sy.nr.tolist()}, {NG} G '
          f'(|q+G|^2 <= {a.ecut_eps} Ry), Q crystal {np.round(Qcr, 4).tolist()} = q {np.round(qcr, 4).tolist()} '
          f'+ G {GQ.tolist()}')
    # memory: the binned matrices, one omega chunk (chi, eps, eps^-1, anti-Hermitian part, ...), the mapped output
    nb = int(np.ceil(a.emax / a.de)) + 2
    wgrid = np.arange(0.0, a.emax, a.de)
    wc = max(1, min(len(wgrid), int(2.5e8 // (16 * NG * NG * 5))))
    need = (16 * (nb + len(wgrid)) * NG * NG + 6.4e7 + wc * 16 * NG * NG * 5
            + (16 * len(wgrid) * NG * NG if a.full else 0)) / 1e9
    print(f'chi0_ref: memory estimate {need:.2f} GB ({nb} bins x {NG}^2, omega chunks of {wc})')
    if need > a.max_mem:
        sys.exit(f'chi0_ref: {need:.2f} GB exceeds --max-mem {a.max_mem} GB; raise --de, lower --emax or --ecut-eps')
    A = binned_A(sy, qcr, Gm, a.emax, a.de)
    herm = max(np.abs(A[b] - A[b].conj().T).max() for b in range(nb))
    mineig = min(np.linalg.eigvalsh(0.5 * (A[b] + A[b].conj().T)).min() for b in range(nb))
    print(f'chi0_ref: max |A - A^H| {herm:.2e}, min eigenvalue of the bins {mineig:.2e}')
    full = None
    if a.full:
        full = np.lib.format.open_memmap(a.full, mode='w+', dtype=complex, shape=(len(wgrid), NG, NG))
        with open(a.full + '.G.txt', 'w') as g:
            g.write('# eps^-1[w, i, j] = eps^-1_{G_i G_j}(q, w), chi_GG\' = sum rho*(G) rho(G\'), '
                    'rho(G) = <c,k+q|e^{i(q+G)r}|v,k>\n')
            g.write(f'# q crystal {qcr.tolist()}  w = 0 .. {a.emax} eV step {a.de}\n# i  G (crystal, b basis)  |q+G|^2 (bohr^-2)\n')
            for i, (gg, qq) in enumerate(zip(Gm, q2)):
                g.write(f'{i} {gg[0]} {gg[1]} {gg[2]} {qq:.8f}\n')
    eQQ = np.empty(len(wgrid), complex); einvQQ = np.empty(len(wgrid), complex)
    ahmin, ahmax = np.inf, 0.0
    chi = chi_kk_fft(sy, A, a.de, len(wgrid), a.eta)    # Kramers-Kronig by FFT on the uniform grid
    del A
    sv = np.sqrt(4 * np.pi / q2)
    for i0 in range(0, len(wgrid), wc):
        ww = wgrid[i0:i0 + wc]
        eps = np.eye(NG)[None] - sv[None, :, None] * chi[i0:i0 + wc] * sv[None, None, :]
        einv = np.linalg.inv(eps)
        eQQ[i0:i0 + wc] = eps[:, iQ, iQ]
        einvQQ[i0:i0 + wc] = einv[:, iQ, iQ]
        # retarded response: the anti-Hermitian part of -eps^-1 is positive semidefinite for w > 0
        ah = -(einv - np.conj(np.swapaxes(einv, 1, 2))) / 2j
        ev = np.linalg.eigvalsh(ah[ww > 0])
        if ev.size:
            ahmin = min(ahmin, ev.min()); ahmax = max(ahmax, ev.max())
        if full is not None:
            full[i0:i0 + wc] = einv
    if full is not None:
        full.flush(); del full
    print(f'chi0_ref: anti-Hermitian part of -eps^-1 over w > 0: eigenvalues in [{ahmin:.2e}, {ahmax:.2e}]')
    L_rpa = -einvQQ.imag
    L_ipa = -(1.0 / eQQ).imag
    im_eps0 = eQQ.imag      # = -v Im chi0_QQ: what turbo_eels calls the IPA loss (no Hartree term)
    epsM = 1.0 / einvQQ
    w = wgrid / HA
    wp2 = 4 * np.pi * sy.nelec / sy.omega
    np.savetxt(a.out, np.column_stack([wgrid, im_eps0, L_ipa, L_rpa, epsM.real, epsM.imag]),
               header='w(eV)  Im eps_NLF(=turbo IPA)  L_NLF  L_RPA  Re eps_M(RPA)  Im eps_M(RPA)')
    print(f'chi0_ref: Im eps without local fields (turbo IPA): peak {wgrid[np.argmax(im_eps0)]:.2f} eV ({im_eps0.max():.3f}), '
          f'f-sum fraction {trapz(w * im_eps0, w) / (0.5 * np.pi * wp2):.3f}')
    for name, L in (('NLF', L_ipa), ('RPA', L_rpa)):
        fs = trapz(w * L, w) / (0.5 * np.pi * wp2)
        print(f'chi0_ref: {name} loss peak {wgrid[np.argmax(L)]:.2f} eV ({L.max():.3f}), f-sum fraction {fs:.3f}')
    print(f'chi0_ref: Re eps_M(w=0) RPA {epsM[0].real:.3f}, no local fields {eQQ[0].real:.3f}')


if __name__ == '__main__':
    main()
