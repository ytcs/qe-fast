# Loss functions and GW: what the tree can do (issue #1)

Issue #1 asks for a route to the impact-ionisation kernel of the ε_eh / Fano-factor calculation: quasiparticle
bands, and the inverse dielectric matrix ε⁻¹_GG′(q, ω) up to about 100 eV. This directory records what the tree
provides today, one measured smoke test, and the recommended route.

## What builds from this tree

Built and linked on CPU (gfortran, CMake, no MPI): `turbo_eels.x` and `turbo_spectrum.x` (TDDFPT) and `pw2bgw.x`
(the BerkeleyGW interface). Neither uses the qe-fast GPU engine. `EPW/examples/lif_hdf5` is an upstream example
that chains QE, BerkeleyGW (epsilon, sigma, kernel, absorption) and EPW for LiF.

| route | gives | limits |
|---|---|---|
| `turbo_eels` (Liouville–Lanczos TDDFPT) | −Im ε⁻¹(Q, Q; ω) for any Q, inside or outside the first zone, with local fields in the induced density (RPA or ALDA) | one diagonal element per Lanczos chain as shipped; NC, US and noncollinear supported upstream |
| BerkeleyGW via `pw2bgw.x` | full ε⁻¹_GG′(q, ω), G0W0 / GW0 quasiparticles, BSE | external code (BSD-3); NC pseudopotentials only in the pw2bgw path; GPU build in BerkeleyGW 4.0 |
| GWW (upstream `pw4gww`, `gww`) | GW with Lanczos polarizability | Γ-point (supercell) formulation, not a route for k-resolved kernels; removed from qe-fast |

The published benchmark the issue cites (Nielsen and Fischetti, arXiv:2610.16949) uses `turboEELS` for the loss
function, with an isotropic loss function in q, and VASP only for the GW/BSE correction of Si.

## Si smoke test

`si-eels/run.sh <bin> <workdir>`: Si, LDA (`Si.pz-vbc.UPF`), 20 Ry, 8×8×8 k, ALDA, 1500 Lanczos steps extrapolated
to 20000, 0.02 eV Lorentzian, spectrum to 150 eV. Q is along x in units of 2π/a.

| Q | loss peak | f-sum violation (0–150 eV) | Re ε(Q, ω = 0) |
|---|---|---|---|
| 0.1 | 17.8 eV | 5.0 % | 12.7 |
| 2.1 (outside the first zone) | 29.9 eV, broad | 6.1 % | 1.65 |

- The small-Q plasmon at 17.8 eV compares with 16.7 eV measured and with about 17.9 eV in the paper above.
- At Q = 0.1 the violation is 5.03 % with both 1000 and 1500 Lanczos steps, so the chain is converged. It is
  4.86 % with a 0.005 eV Lorentzian and 4.05 % when the integral runs to 400 eV. What makes up the remainder has
  not been separated; the nonlocal pseudopotential, for which the f-sum rule is not exact, is one candidate.
- Q = (2.1, 0, 0) and (0, 2.1, 0) agree to 1.5×10⁻⁹. That checks the symmetry handling only.
- The table is the output of `run.sh` as committed. One Q took 4–5 minutes on 7–14 OpenMP threads of a desktop CPU.

## Inverse dielectric matrix by sum over bands

The impact-ionization rate of the book's E-13 contracts the full Im W_GG′(q, ω) with pair densities, so the
off-diagonal elements are needed, at every energy transfer up to the carrier energy. The route is the sum over
bands that E-16 specifies (the independent-particle matrix R⁰_GG′, then a Dyson solve), which gives the whole
matrix at every ω in one pass. A Lanczos column per G (extending `turbo_eels`) would cost on the order of
10⁸–10⁹ H applications per q for Si against about 10⁵ pair-density FFTs, and stays as an independent check of
the diagonal.

`qefast/gw/chi0_ref.py` is the CPU reference: pair densities ρ(G) = ⟨n′, k+q| e^{i(q+G)r} |n, k⟩ for any band pair,
transitions binned in energy (each bin's matrix is Hermitian and positive semidefinite; the run checks both),
χ⁰(ω) by a Kramers–Kronig sum, ε = 1 − v^½ χ⁰ v^½, and ε⁻¹ by inversion. Norm-conserving pseudopotentials only.
`si-chi0/run.sh` runs it on the Si of the smoke test (8³ k without symmetry, 60 bands, local-field basis
|q+G|² ≤ 8 Ry, 0.1 eV Lorentzian); `APPROX=IPA` and `APPROX=RPA_with_CLFE` of `si-eels/run.sh` give the
`turbo_eels` side at the same Q, which must be on the k grid (Q = 0.25 and 2.25 along x are; 0.1 and 2.1 are not).

Static values (8 Ry basis, 60 bands):

| Q (2π/a, x) | Re ε(Q, 0) without local fields | with local fields |
|---|---|---|
| 0.25 | 10.51 (turbo_eels 10.48) | 9.48 (9.46) |
| 2.25 | 1.459 | 1.408 (1.410) |

Spectra, 120 bands, both codes at turbo_eels' 0.02 eV Lorentzian, then the same Gaussian smoothing (σ = 1 eV) of
both (an 8³ grid leaves spikes at 0.02 eV; turbo_spectrum's own larger broadenings change the integral and are not
used):

| Q | quantity | peak, chi0_ref / turbo_eels | integral 0–58 eV, ratio | largest difference / peak |
|---|---|---|---|---|
| 0.25 | Im ε without local fields (turbo's IPA "loss") | 4.42 / 4.44 eV | 1.039 | 0.17 |
| 0.25 | RPA loss, local fields | 17.20 / 17.14 eV | 1.014 | 0.09 |
| 2.25 | Im ε without local fields | 23.88 / 23.78 eV | 1.014 | 0.03 |
| 2.25 | RPA loss, local fields | 31.78 / 31.78 eV | 1.028 | 0.05 |

- Bands: at Q = 0.25 the loss peak without local fields is 17.08, 17.06 and 17.06 eV with 20, 40 and 60 bands; at
  Q = 2.25 the RPA f-sum fraction is 0.54, 0.86 and 0.93 with 30, 60 and 120 bands, and the 0.1 eV peak stays at
  32.72 eV from 60 to 120 bands. The f-sum fraction at Q = 0.25 is 0.88, 0.90, 0.91 and 0.94 with 20, 40, 60 and
  120 bands, the same with and without local fields; `turbo_eels` misses it by 5.6 %.
- The first version conjugated the antiresonant term's matrix. The head is real, so the diagonal checks passed and
  the local-field effect came out half as large (Re ε(0) 9.98 against 9.46); with time reversal both terms carry the
  same matrix.
- Memory: the binned matrices are 16 × bins × N_G² bytes (0.6 GB for 3000 bins and 109 G); the script estimates its
  use first and stops above `--max-mem` (2 GB by default).

### chi0.x (host path)

`chi0.x` (`qefast/gw/chi0.f90`, built with the PP programs) is the same calculation inside the QE build: it reads the
run with QE's own modules, uses QE's FFTs, one ZGEMM per energy bin and LAPACK for the inversion, sizes its energy
chunks from a memory budget (`max_mem`) and writes the same columns, plus the whole ε⁻¹ matrix on request
(`filfull`, with the G list). It runs on one MPI process for now. Against `chi0_ref.py`:

- Si, 8³ k, 60 bands, 8 Ry basis, Q = 0.25: the whole ε⁻¹ (400 energies × 109²) agrees to 2×10⁻¹³; the output
  columns agree to their printed precision.
- CI runs `si-chi0/ci.sh`: Si on 4³ k with 12 bands at Q = 0.5, whole matrix and columns to 10⁻⁸ (measured
  1.9×10⁻¹² and 5×10⁻¹⁰). Its first run caught a real error: QE's `'Wave'` transform covers only the columns of
  the wavefunction sphere, and the k+q coefficients shifted by an umklapp vector fell partly outside them; the
  program now uses the full-box transform.
- Kramers–Kronig is a convolution and a correlation on the uniform energy grid, done by FFT in both codes (equal to
  the direct sum to 2×10⁻¹⁵). On the Si 8³ case (one thread) it cut that stage from 84 to 19 CPU s and the run from
  150 to 89 CPU s; pair densities (38 s) are now the largest stage.

### chi0.x on the GPU

In CUDA builds `chi0.x` runs the same stages on the device (`qefast/gw/chi0_dev.f90`): pair densities for all
conduction bands of a (k, v) in one batched cuFFT, transitions buffered on the device and binned with one cuBLAS
ZGEMM per occupied bin, Kramers–Kronig by strided 1D cuFFTs, and batched cuBLAS LU for ε⁻¹. It falls back to the host
path when one k-point does not fit in device memory (`device = .false.` forces the host path). H100, Si 8³ k, 60
bands, Q = 0.25:

| run | whole ε⁻¹ against chi0_ref.py | program wall time (s) |
|---|---|---|
| device, 8 Ry (109 G), 400 energies | 1.7×10⁻¹³ | 5.5 |
| host path in the same GPU build (8 OpenMP threads) | 2.0×10⁻¹³ | 16.3 |
| device, 8 Ry, 3000 energies to 60 eV | — | 5.4 (reading the wavefunctions 4.0) |
| device, 16 Ry, 3000 energies | — | 8.1 |

The same 8 Ry, 3000-energy case takes 89 CPU s on one desktop core with the host path. On the device, reading and
transforming the wavefunctions on the host is now the largest stage. The PSD check of −ε⁻¹ runs on 100 energies by
default (`check_psd`); on every energy of the 16 Ry case it took 30 s on the host.

Two errors found on the way, both in the host path of GPU builds: QE's internal FFT backend (used there) scales a
batched 1D transform differently from a single one and transforms its input in place, so the program measures the
scalings per batch size on a fresh impulse.

## Impact-ionization rate (E-13): contraction reference

`qefast/gw/gamma_ref.py` contracts Im W with pair densities, Im Σ_nk = (1/(N_k Ω)) Σ_{q≠0, n′} φ† Im W(q, ε_nk − ε_n′k−q) φ
with φ(G) = ⟨n,k| e^{i(q+G)r} |n′, k−q⟩ in the orientation of χ⁰, and Γ_ii = 2z|Im Σ|/ħ. It shows what the program
has to produce: the full ε⁻¹_GG′(q, ω) is never stored. For Si at the E-13 settings (about 29 q in the irreducible
wedge and 512 in the zone, 1000 energies, 200 G) it would take 19 GB on the wedge and about 330 GB on the zone; the
program contracts per q and keeps Γ_ii and the final-state table.

Si, 4³ k without symmetry, 30 bands, 4 Ry basis, W with a 0.05 eV Lorentzian, z = 1, q = 0 term left out (1/64
of the zone); LDA gap on this grid 0.64 eV:

| E − E_cbm (eV) | 0–0.5 | 0.5–1 | 1–1.5 | 2.5–3 | 4–4.5 | 5.5–6 |
|---|---|---|---|---|---|---|
| mean Γ_ii (s⁻¹) | 1.6×10¹¹ | 1.3×10¹² | 1.5×10¹³ | 1.3×10¹⁴ | 4.3×10¹⁴ | 7.0×10¹⁴ |

- Im Σ ≤ 0 for every state, as it must be for a retarded W.
- Below the gap energy the rate must vanish. The 1.6×10¹¹ s⁻¹ there is the Lorentzian tail of Im W below the gap:
  with a 0.01 eV Lorentzian it falls to 3.0×10¹⁰ s⁻¹, in proportion to the width.
- This is a prototype grid: the conduction minimum near X is not on 4³, the LDA gap is half the measured one, and
  the q = 0 term is missing. It fixes the output and the orientation, not the numbers.

### gamma_ii.x (host path)

`gamma_ii.x` computes the same rates: per q it bins the transitions, builds ε⁻¹ by Kramers–Kronig and inversion
(the `chi0.x` stages, shared through `qefast/gw/gw_common.f90`), turns it into Im W on the energy grid, and
contracts Im W with the pair densities of the hot carriers. Only one q is held in memory at a time. Input
`&inputgii`: `prefix, outdir, emax_carrier` (eV above E_cbm), `ecut_eps` (Ry), `de, eta` (eV), `scissor` (eV),
`fileout` (band, k, E − E_cbm, Im Σ, Γ_ii), `filpart` (the rate into every final state (n′, k′) with both energies;
the rates into the final states of a carrier sum to its Γ_ii), `max_mem` (GB).

Si, 4³ k without symmetry, 30 bands, carriers up to 2 eV, 4 Ry basis, 0.02 eV grid, 0.05 eV Lorentzian, two
threads on a shared desktop CPU:

| | Γ_ii per state | partial rates | time |
|---|---|---|---|
| `gamma_ref.py` | – | 1929 | 350 s |
| `gamma_ii.x` | equal to 7 digits (57 of 63 states; the other 6 are tail values below 10⁻⁶ of the peak) | same 1929 keys, equal to 8×10⁻¹⁰ of the peak | 188 s |
| `gamma_ii.x`, H100 | equal to 9×10⁻⁵ (from its own nscf on the GPU) | – | 3.8 s |

Three quarters of the host time is the pair-density FFTs and a fifth the wavefunction reads, both repeated for every
q. In CUDA builds `gamma_ii.x` reads the plane-wave coefficients of every k once, keeps them on the device, and runs
every stage there: periodic parts and pair densities (batched cuFFT), binning (cuBLAS ZGEMM per bin),
Kramers–Kronig, the batched inversion to Im W, and the contraction (one ZGEMM per grid energy over the pairs that
interpolate at it). `device = .false.` keeps the host path. On the H100 the device path equals the reference on the
CI case (12 bands; Γ_ii to 7 digits, 1916 partial rates to 8×10⁻¹⁰ of the peak) and on a noncollinear Si 4³ run
(24 spinor bands; 126 states, 7664 partial rates to 10⁻⁹ of the peak) in 4.4 s, against 114 s for the host path on
8 cores of the same machine. CI runs the comparison on 12 bands (carriers to 2.5 eV) in `si-chi0/ci.sh`.

Stages from here to the issue's gate, each with its own check:

1. ε⁻¹ per q: the reference above, then the GPU program element-wise against it (full matrix, PSD check).
2. Γ_ii by contraction, including the q = 0 head and wings: positivity, threshold (zero below E_g), convergence in
   k grid, bands and basis; then a comparison with published GW impact-ionization rates of Si.
3. Phonon emission and absorption rates from EPW for the same states.
4. The cascade Monte Carlo of chapter 18 with both kernels: ε_eh and F_eh at 123 K against 3.749 eV and 0.119.

## Recommended route

1. **ε⁻¹_GG′(q, ω): native, sum over bands.** The reference above, then a GPU program checked against it.
2. **Quasiparticle corrections: BerkeleyGW** through `pw2bgw.x`, for the G0W0 shifts and lifetimes. These are
   needed once per material, not per Q.
3. **Excitons in LiF and NaF: BerkeleyGW BSE.** Yambo is GPL and would also do; neither has been adopted, and
   whichever is adopted goes into NOTICE first.

## Open questions

- The issue's Si gate, ε_eh = 3.749 ± 0.001 eV and F_eh = 0.119 ± 0.002, is a measurement at 123 K (Rodrigues et
  al. 2021). The paper above compares with 3.6 eV, a room-temperature value, and gives no Fano factor. A gate
  must state its temperature.
- The 100 eV range reaches the Si 2p edge (about 99 eV), which `Si.pz-vbc.UPF` does not have. A gate up to that
  range needs a pseudopotential with the 2s2p shell in the valence.
- The issue's targets are insulators, and CaWO₄ needs spin–orbit coupling. The batched ph.x engine declines both
  (smearing and collinear only), so the EPW step of the pipeline runs QE's own loop for them today.
