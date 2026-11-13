# qe-fast switches

Every run-time switch of qe-fast is an environment variable read through `qef_getenv`
(`qefast/config/qef_config.f90`). Its value is, in order:

1. the environment variable, when it is set;
2. the entry of the active preset in the table of `qef_config`;
3. unset: the code's own default for that switch.

`QEF_PRESET` is `off`, `safe`, `fast` or `max`. Any other value stops the run. When it is unset, the
built-in default applies: **`safe`**. Each program prints the active preset and any table override at
start (`qe-fast preset:`, `qe-fast override:`). `mpirun` must pass the variables on to remote ranks;
`qefast/examples/eph-tc/bench/run-chain.sh` runs a whole chain under one preset and passes every `QEF_*`.

## Classes

- **B**: bitwise identical output.
- **E**: same mathematics in another order: differences at rounding level, inside the gate.
- **N**: changes the numbers: the answer moves, and the gate decides whether that is acceptable.

The gate, against the same material's `off` run: every mode at every q within max(0.5 %, 1 cm⁻¹), the three
acoustic modes at Γ excluded; λ and ω_log within 1 %; Allen–Dynes T_c within max(2 %, 0.5 K). B1 is MgB₂ and
B2 is Mg₂IrH₆, the two tuning benchmarks. G1–G8 are the eight held-out materials: MgO, Al, Nb, H₃S at 200 GPa,
LiAl₂Rh, TiN (39³ FFT grid forced), MgB₃C₃, B₁₃C₂. All gate runs are GPU runs (A40, H100); no preset has been
gated on a CPU build.

## Presets

| preset | what it is | status |
|---|---|---|
| `off` | The engine with every lever off (`QEF_STH=1 QEF_STRICT=0 QEF_STH_DVB=0 QEF_W90_NT=1`). Not plain QE: `QEF_STH=0` gives QE's own loop. | Reference for every gate. Bitwise identical to the earlier all-off path: 21 files on B1 smoke (dyn, dvscf, a2f) against an earlier binary; 34 files on G2, G3, G5, G6 against the archived all-off runs. |
| `safe` | All levers of classes B and E that leave the Wannier gauge alone. | **Full gate passed** on G2, G3, G5, G6, the four materials that failed the study's earlier `safe`, and on B1 and B2 smoke. Not run itself on G1, G4, G7, G8: they passed the frequency gate with that earlier, larger set, and G1, G4, G7 its full gate; G8 has no λ reference. |
| `fast` | The held-out study's `fast` set: `safe` + `QEF_PW_SYMROT=1 QEF_W90_WNT=8 QEF_W90_MMN=1` + lever C (`QEF_THRESH_L=10 QEF_STH_DEFL=1`). | **Phonon-gated only on held-out materials.** Frequencies pass on G1–G8. λ, ω_log or T_c fail on G2, G3, G5, G6; G8 has no λ reference. |
| `max` | The held-out study's `max` set: `safe` + `QEF_PW_SYMROT=1 QEF_W90_WNT=8 QEF_W90_MMN=1` + CBSR (`QEF_STH_CBSR=1 QEF_STH_CBSR_XCH=1 QEF_STH_CBSR_HOST=1`, automatic check on). | **Phonon-gated only on held-out materials.** Same failures as `fast`. |

`fast` and `max` keep the definitions they were gated with, so they carry the three switches `safe` leaves out,
and their λ is subject to the gauge hazards below. To drop those from either preset, set
`QEF_PW_SYMROT=0 QEF_W90_WNT=1 QEF_W90_MMN=0`; that combination has not been gated.

Measured speedup, end to end against `off`, median / worst over the seven held-out materials with a total
(A40): gen `safe` 1.38× / 1.11×, `fast` 1.68× / 1.42×, `max` 1.95× / 1.09×. Those `safe` figures include the three
switches now left out. The default `safe` is slower than that, because the three switches carried the nscf and part of the Wannier gain:

| | G2 Al | G3 Nb | G5 LiAl₂Rh | G6 TiN | B1 smoke | B2 smoke |
|---|---|---|---|---|---|---|
| `safe` end to end | 1.50× | 1.35× | 1.07× | 1.10× | 1.50× | 1.61× |
| largest \|Δω\| (cm⁻¹) | 0.054 | 0.150 | 0.000 | 0.002 | 0.016 | 0.011 |
| rel Δλ | −6.3e-5 | −4.8e-4 | −1.6e-6 | −9.1e-6 | 3.0e-6 | 7.7e-6 |
| largest Cartesian WF-centre shift (Å) | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |

Margin used on G2, G3, G5, G6: 3.8 %, 15.0 %, 0.0 %, 0.2 %.

## Things to know before you rely on this

- **FFTDX and FDXDP are compiled for two grids only**: 30×30×36 and 54×54×54 (`FDX_GRIDS` in
  `qefast/sth/fdx/sth_fftdx.cu`). They declined on 8 of 8 held-out materials. The decline is clean
  (`qef fftdx: unsupported grid/box`) and the run keeps cuFFT, but the largest lever in `safe` then does nothing:
  `safe` left the solve stage where it was on five of seven held-out materials.
- **`QEF_STH_DVB=1` has no fallback once the engine declines at run time.** There is no dV_bare buffer to go back
  to, so a decline, for example for lack of device memory, stops with `engine declined with QEF_STH_DVB=1`, even
  with `QEF_STRICT=0`. `safe`, `fast` and `max` all set it. Set `QEF_STH_DVB=0` for a run that stops this way.
- **`QEF_EPW_NOVME=1` stops if the fused fine path declines at run time** (pair-kernel shared memory, images)
  in an `ephwrite` + Eliashberg run: the velocity matrices were never built. Set `QEF_EPW_NOVME=0` for such a
  run. The switch turns itself off in CPU builds and when an input flag rules out the fused path or reads
  velocities; before 2026-11-12 it stopped those runs too, including every CPU Eliashberg run under `safe`.
- **B1's Wannierisation ends on a mirror saddle.** At the recipe's endpoint every WF centre sits on the σh plane
  and the spread has a negative curvature (θ_min −0.0021, ‖g‖ 6.5e-8), so wannier90 stops converged on a strict
  symmetry-protected saddle (Ω 7.617 Å² on the frozen chain; the minimum is Ω 7.581 Å² on the proxy chain, with
  λ moved by 6.3e-4 relative). B2 has its own saddle at Ω 41.935 Å² until about iteration 370. Whether a run
  leaves such a point depends on rounding, so a rounding-level change can move λ.
- **λ depends on the home cells of the WFs at unchanged spread.** Moving individual WFs by lattice vectors is a
  pure gauge change: the spread does not move, and λ moved by +3.7 % on the B2 nq 2³ proxy (the size at nq 4³ was
  not measured). Identity checks must compare **Cartesian centres**, not spread sums
  (`qefast/examples/eph-tc/bench/gen/wfcentres.py`).
- **The `cmp`-identical `max` against `safe` on G1, G7 and G8 proves the fallback path only.** These are the three
  materials where the CBSR check declined on payoff, so `max` ran the CG path there. It says nothing about CBSR
  itself, which is class N and gated on frequencies.
- **Lever C was tuned on B1 and B2.** Its worst held-out frequency margin is 38.7 % (G7), against 60 % used on B2.
  It has passed the frequency gate on all eight held-out materials and no held-out λ gate.
- **Engine scope.** Magnetic, SOC, USPP/PAW and DFT+U systems, `occupations = 'fixed'` insulators and runs with more
  than one band group take QE's own loop, with no speedup.
- **SOC-W (certified Wannier minimisation) is not included.** v2 failed on B2.

## Engine

| switch | default | class | gate evidence | fallback |
|---|---|---|---|---|
| `QEF_STH` | 1 in all presets | — | batched Sternheimer engine, batched nscf and GPU tail; the baseline of the H100 benchmarks | 0: QE's own loop. Unsupported cases take it automatically. |
| `QEF_STRICT` | 0 | — | 1 turns any decline into a stop; used for timing runs | — |
| `QEF_TAIL`, `QEF_NSCF_BATCH`, `QEF_NSCF_ROT` | on (0 disables) | E | part of `off` | 0: QE's own code for that step |

## ph.x levers (on in `safe`, `fast`, `max`)

| switch | safe | class | gate evidence | fallback |
|---|---|---|---|---|
| `QEF_STH_DVB` | 1 | E | B1, B2 smoke chains PASS; frequencies pass G1–G8 | **none** at run time, see above |
| `QEF_STH_FFTDX` | 1 | E | B1 0.081, B2 0.000 cm⁻¹; H100 ph.x 1.47× / 1.62× | cuFFT on any other grid (8 of 8 held-out) |
| `QEF_STH_FDXDP` | 1 | E | B1 frozen q2 0.0077 cm⁻¹ with OBATCH | needs FFTDX; else cuFFT |
| `QEF_STH_OBATCH` | 1 | E | B1 frozen q2 0.0077 cm⁻¹ with FDXDP; q7 dyn and dvscf `cmp`-identical under QCACHE | — |
| `QEF_STH_QCACHE` | 1 | B | B1 frozen q7 `cmp`-identical | — |
| `QEF_NSCF_SYMROT` | 1 | E | B1, B2 gate chains PASS (≤ 0.2 cm⁻¹, λ ≤ 0.2 %); shared-gauge run clean on G3, G5, G6 | declines for noncollinear, Γ-only, LSDA, DFT+U, USPP, more than one pool or band group; unconverged representative: its orbit is diagonalised |
| `QEF_PHQ_FUSE` | 1 | E | with NSCF_SYMROT, same gate | — |
| `QEF_PHQ_VKB` | 1 | B | bitwise against SYMROT + FUSE; needs FUSE | — |
| `QEF_DVXC_GPU` | 1 | E | B1, B2 smoke ≤ 0.2 cm⁻¹, λ ≤ 0.2 % | CPU path |
| `QEF_STH_KDEV` | 2 | B | bitwise, B1 frozen q7 | 0: host loop |
| `QEF_STH_PSIR_PART` | 1 | E | dvscf 5e-14; B2 smoke chain PASS | — |
| `QEF_STH_DRHO2` | 1 | E | dvscf rel L2 1.5e-15, frequencies identical | — |

Every lever above is also covered by the shared-gauge diagnostic: with all-off's Wannier output, each set's ph.x data
moved λ by at most 4.7e-4 relative on G3, G5, G6.

## wannier90 and EPW levers (on in `safe`, `fast`, `max`)

| switch | safe | class | gate evidence | fallback |
|---|---|---|---|---|
| `QEF_W90_NT` | 8 | B | `dis_extract` over k; G6 λ unchanged to all digits; 1 in `off` | 1: serial |
| `QEF_W90_AMN` | 1 | B | B1, B2 gate; G6 λ unchanged to all digits | — |
| `QEF_W90_NOHR` | 1 | B | every gated quantity identical | — |
| `QEF_EPW_FDX` | 1 | E | full gate B1, B2; EPW group leaves G6 λ unchanged; `_CHECK` prints the difference | cuFFT path off grid |
| `QEF_EPW_BEC` | 1 | E | as above | — |
| `QEF_EPW_STAR` | 1 | E | as above; `QEF_EPW_STAR_FB=n` tunes its fallback | direct computation of deferred k |
| `QEF_EPW_B2W` | 1 | E | epmatwp rel L2 3e-16 to 5e-16, λ identical | — |
| `QEF_EPW_IDLEWAIT` | 1 | B | `.ukk` and epmatwp identical | — |
| `QEF_EPW_NOVME` | 1 | B | `.ukk` and epmatwp `cmp`-identical | **none**, see above |

## Opt-in: in `fast` and `max` as gated, never in `safe`

| switch | fast, max | class | evidence | why not in safe |
|---|---|---|---|---|
| `QEF_PW_SYMROT` | 1 | N | G6: spread 11.93 → 11.67 Å², λ 0.906 → 1.173 | Reaches a different Wannier minimum. Its λ is not wrong by itself (all-off is not a truth reference for λ), but it is not the same answer. |
| `QEF_W90_WNT` | 8 | N | G6 alone: λ −1.8 % at unchanged spread (0.8895) | Not bitwise (trajectory differs at 1e-10), and the saddles above turn that into λ. |
| `QEF_W90_MMN` | 1 | N | G6 alone: λ −1.8 % at unchanged spread (0.8895) | As WNT. |

## Other opt-in levers

| switch | presets | class | evidence | note |
|---|---|---|---|---|
| `QEF_THRESH_L=10`, `QEF_STH_DEFL=1` (lever C) | `fast` | N | ph.x 1.21× / 1.23× on B1 / B2; B2 margin used 60 %; held-out 38.7 % | tuned on B1 and B2 |
| `QEF_STH_CBSR=1` | `max` | N | B1 smoke 0.126 cm⁻¹; held-out frequencies pass | exact spectral solve. `QEF_STH_CBSR_AUTO` is on whenever CBSR is set: it checks memory and payoff once at the first q ≠ Γ and falls back to CG. Under `QEF_STRICT=1` it only reports. |
| `QEF_STH_CBSR_XCH`, `QEF_STH_CBSR_HOST` | `max` | B | with CBSR | memory path of the check |
| `QEF_STH_CBSR_RF` and `_TF32`, `_ACT`, `_M0`, `_CMAX`, `_RFMAX`, `_CHK`, `_PAD`, `_GAMMA` | none | N | refinement form; B2 margin used 21 / 56 / 126 % | research only |
| `QEF_EPW_NOWRITE=1` | none | B | λ identical | skips epmatwp and `.epb`, which breaks later EPW-only reuse (screening, shared-gauge runs) |
| `QEF_STH_MP=1`, `QEF_MP_ETA`, `QEF_MP_SAFE` | none | N | B1, B2 PASS, 74 % margin used on B2 with lever C; about 2× on ph.x on a 4090 | FP32 inner CG; hardware dependent |
| `QEF_TR2_SCALE` | none | N | 0.01 passes against a tight same-grid reference at B2 q2, q3 with exact CBSR | tightens tr2_ph |
| `QEF_STH_DRF=1`, `QEF_EPW_B2WG=1` | none | E | dvscf 1.6e-15 (DRF); no A40 gain | H100 A/B not run |
| `QEF_STH_ETA` | none | N | not adopted | leave unset; wired through the engine's CG kernels, to go with the next engine change |
| `QEF_STH_NOPSIR=1` | none | B | B2 q5 engine memory 68.4 → 24.1 GiB, +10 % solve | memory relief |

## Diagnostics: all opt-in, in no preset

`QEF_STH_THRESH_SCALE` and `QEF_STH_TRUERES` change the answer. `QEF_STH_MEMLOG`, `QEF_STH_TIME`, `QEF_STH_PROF`,
`QEF_EPW_PROF`, `QEF_STH_DUMP`, `QEF_STH_DUMP_CALL` and the `*_CHECK` switches (`QEF_STH_FFTDX_CHECK`,
`QEF_STH_FDXDP_CHECK`, `QEF_STH_DVB_CHECK`, `QEF_NSCF_SYMROT_CHECK`, `QEF_PW_SYMROT_CHECK`, `QEF_W90_MMN_CHECK`,
`QEF_EPW_FDX_CHECK`, `QEF_EPW_STAR_CHECK`, `QEF_EPW_ELIASH_CHECK`) only print. `QEF_NSCF_SYMROT_TIGHT`,
`QEF_PW_SYMROT_TIGHT`, `QEF_NSCF_SYMROT_PREFER_FT`, `QEF_NSCF_SYMROT_REVOPS`, `QEF_STH_QCACHE_NOKQG`, `QEF_STH_PSIR_RES`,
`QEF_STH_PSIR_KC`, `QEF_STH_CBSR_MINGAIN`, `QEF_STH_CBSR_RESERVE` and `QEF_STH_CBSR_HOSTN` are test and tuning knobs.
`QEF_W90_DUMP` exists only as a saved patch of the SOC-W study, not in this tree.

## Default-on EPW paths (part of `off`; 0 restores QE's code)

`QEF_EPW_COARSE`, `QEF_EPW_NL`, `QEF_EPW_FINE`, `QEF_EPW_ELIASH`, `QEF_EPW_ELIASH_KERNEL`, `QEF_EPW_ELIASH_MIX`,
`QEF_EPW_ELIASH_WARM`, `QEF_EPW_QBATCH`. `QEF_EPW_VME=1` forces the velocity matrices; `QEF_EPW_ETFMEM=keep` keeps the
input's `etf_mem`. All were gated in the full H100 baseline chain.

## Script and build switches

`QEF_RUN_DIRECT`, `QEF_RUN_PREWARM`, `QEF_RUN_NOEPW` (read by `run-chain.sh`, with `QEF_EPW_GOFILE` for the warm start).
CMake: `QEF_MATHDX_INC`, `QEF_CUTLASS_INC` (without both, FFTDX is a stub that always declines), `QEF_NVCC`.
`FDX_*` variables are timing knobs of the FFTDX kernels.
