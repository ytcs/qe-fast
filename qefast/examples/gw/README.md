# Loss functions and GW: what the tree can do (issue #1)

Issue #1 asks for a route to the impact-ionisation kernel of the ε_eh / Fano-factor calculation: quasiparticle
bands, and the inverse dielectric matrix ε⁻¹_GG′(q, ω) up to about 100 eV. This directory records what the tree
provides today, one measured smoke test, and the recommended route.

## What builds from this tree

Built and linked on CPU (gfortran, CMake, no MPI): `turbo_eels.x` and `turbo_spectrum.x` (TDDFPT), `pw2bgw.x` (the
BerkeleyGW interface), and `pw4gww.x` and `gww.x` (GWW). None of them uses the qe-fast GPU engine.

| route | gives | limits |
|---|---|---|
| `turbo_eels` (Liouville–Lanczos TDDFPT) | −Im ε⁻¹(Q, Q; ω) for any Q, inside or outside the first zone, with local fields in the induced density (RPA or ALDA) | one diagonal element per Lanczos chain as shipped; NC, US and noncollinear supported upstream |
| BerkeleyGW via `pw2bgw.x` | full ε⁻¹_GG′(q, ω), G0W0 / GW0 quasiparticles, BSE | external code (BSD-3); NC pseudopotentials only in the pw2bgw path; GPU build in BerkeleyGW 4.0 |
| GWW (`pw4gww`, `gww`) | GW with Lanczos polarizability | Γ-point (supercell) formulation; not a practical route for k-resolved kernels |

The published benchmark the issue cites (Nielsen and Fischetti, arXiv:2610.16949) uses `turboEELS` for the loss
function, with an isotropic loss function in q, and VASP only for the GW/BSE correction of Si.

## Si smoke test

`si-eels/run.sh <bin> <workdir>`: Si, LDA (`Si.pz-vbc.UPF`), 20 Ry, 8×8×8 k, ALDA, 1500 Lanczos steps extrapolated
to 20000, 0.02 eV Lorentzian, spectrum to 150 eV. Q is along x in units of 2π/a.

| Q | loss peak | f-sum violation (0–150 eV) | Re ε(Q, ω = 0) |
|---|---|---|---|
| 0.1 | 17.8 eV | 5.0 % | 12.7 |
| 2.1 (outside the first zone) | 29.8 eV, broad | 6.1 % | 1.65 |

- The small-Q plasmon at 17.8 eV compares with 16.7 eV measured and with about 17.9 eV in the paper above.
- The violation is the same with 1000 and 1500 Lanczos steps, so the chain is converged. It falls to 4.1 % when
  the integral runs to 400 eV: most of what remains is the Lorentzian tail; the rest is consistent with the
  nonlocal pseudopotential, for which the f-sum rule is not exact.
- Q = (2.1, 0, 0) and (0, 2.1, 0) agree to 1.5×10⁻⁹. That checks the symmetry handling only.
- One Q took 4–5 minutes on 7–12 OpenMP threads of a desktop CPU.

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
- The issue's Si gate (ε_eh = 3.75 eV, F_eh = 0.119) needs a source. The paper above quotes 3.6 eV measured and
  gives no Fano factor.
- The 100 eV range reaches the Si 2p edge (about 99 eV), which `Si.pz-vbc.UPF` does not have. A gate up to that
  range needs a pseudopotential with the 2s2p shell in the valence.
- The issue's targets are insulators, and CaWO₄ needs spin–orbit coupling. The batched ph.x engine declines both
  (smearing and collinear only), so the EPW step of the pipeline runs QE's own loop for them today.
