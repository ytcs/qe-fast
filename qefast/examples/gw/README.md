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

## Recommended route

1. **Loss function: native.** `turbo_eels` already gives the diagonal loss function at arbitrary Q. The EELS
   Lanczos chain is projected onto one vector (ζ); the optical case projects the same chain onto three. Building
   the start vectors exp(i(q+G′)·r) for a set of G′ and projecting onto all of them gives a whole column
   ε⁻¹_G′G(q, ω) per chain at a cost per step of one dot product per G′. That is the extension to make if the
   kernel needs off-diagonal elements. It is also the part worth putting on the GPU.
2. **Quasiparticle corrections: BerkeleyGW** through `pw2bgw.x`, for the G0W0 shifts and lifetimes. These are
   needed once per material, not per Q.
3. **Excitons in LiF and NaF: BerkeleyGW BSE.** Yambo is GPL and would also do; neither has been adopted, and
   whichever is adopted goes into NOTICE first.

## Open questions

- Is the diagonal ε⁻¹(q+G, q+G; ω) enough for the kernel, or are off-diagonal elements required? The answer
  decides whether step 1's extension is needed.
- The issue's Si gate, ε_eh = 3.749 ± 0.001 eV and F_eh = 0.119 ± 0.002, is a measurement at 123 K (Rodrigues et
  al. 2021). The paper above compares with 3.6 eV, a room-temperature value, and gives no Fano factor. A gate
  must state its temperature.
- The 100 eV range reaches the Si 2p edge (about 99 eV), which `Si.pz-vbc.UPF` does not have. A gate up to that
  range needs a pseudopotential with the 2s2p shell in the valence.
- The issue's targets are insulators, and CaWO₄ needs spin–orbit coupling. The batched ph.x engine declines both
  (smearing and collinear only), so the EPW step of the pipeline runs QE's own loop for them today.
