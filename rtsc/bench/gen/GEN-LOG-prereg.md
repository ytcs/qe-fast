# GEN-LOG: held-out generalisation check of the speed levers

Question (user, 2026-10-31): do the speedups hold beyond the two tuned benchmarks (B1 MgB2, B2 Mg2IrH6), or did we overfit?
Agent "gen". Pod `rtsc-gen` (3b8eoeejpdby91): A40 48 GB (sm_86), SECURE EU-SE-1, $0.59/h, 9 vCPU (cgroup quota 7.65 CPUs), 503 GB host RAM, created by the
coordinator at 02:10 UTC 2026-11-06; ssh up 02:14; `nvidia-smi` idle (0 %, 0 MiB, 30 W, no compute apps), driver 595.91.07. Budget ≤ $6.
Fork branch `gen` off origin/rtsc 17883ebc0. Restored /cache/rtsc-b2x-sm86.tar.zst, deps.sh, no stale LR_Modules sth_nscf_symrot.mod.

## Pre-registration (written 2026-11-06 ≈ 02:30 UTC, before any held-out run; never edited after a number)

Nothing below has been run on any held-out material. Runs so far: the Part 1 fallback tests, on B1 and B2 smoke only.

### What was already known when this was written (from reading the code, no run)
- **FFTDX and FDXDP are compiled for exactly two grids.** `FDX_GRIDS` in rtsc/sth/fdx/sth_fftdx.cu is 30×30×36 (B1) and 54×54×54 (B2). Any other grid declines at `fdx_setup` (code 2) and
  keeps cuFFT. So on held-out materials these two levers are expected to do nothing, unless a grid happens to match. This is an overfit by construction. I do not add grids here:
  the question is how the levers as built behave.
- **The engine requires smearing** (`sth_case_ok`: `.NOT. lgauss` returns). A true `occupations = 'fixed'` insulator therefore takes the vanilla loop.
- **QE never picks an FFT size with 7, 11 or 13** in this build (`allowed` in fft_support.f90: only 2, 3, 5). Such a grid exists only if the input sets nr1, nr2, nr3, which bypass `good_fft_order`
  (fft_types.f90). So the awkward-grid material forces it from the input, and that is a stress test, not a case a default run reaches.
- **DVB has no vanilla fallback once the engine declines at run time** (`sth_batched_kernel` stops with 'engine declined with RTSC_STH_DVB=1'). A decline for lack of device memory is therefore a
  crash under `safe`, even with RTSC_STRICT off. The largest cell may or may not reach it.

### Materials (never used for lever tuning)
Screening grade: SG15, the funnel's G5 ecut rule, MP smearing 0.02 Ry, conv_thr 1e-10, tr2_ph 1e-14, the funnel's EPW α²F-only form (`eliashberg = .true.`, no aniso solve, μ* 0.10), Wannier = funnel setup A
(atomic projections, G1 semicore exclusion, frozen window E_F + largest w of {2, 1.5, 1, 0.6} eV with ≤ nWF bands). Grids are trimmed from the funnel's where shown, so each chain costs minutes.

| id | material | atoms, ops | stresses | ecut | scf k | q | EPW k / fine k, q |
|---|---|---|---|---|---|---|---|
| G1 | MgO, rock salt, a 4.25 Å, smearing input | 2, 48 | light-element insulator; frequencies only; no EPW stage (no Fermi surface) | rule | 6³ | 2³ | — |
| G2 | Al, fcc, a 4.04 Å | 1, 48 | 1-atom metal: small N_pw, small batches; launch overheads dominate | rule | 12³ | 4³ | 8³ / 16³, 8³ |
| G3 | Nb, bcc, a 3.31 Å | 1, 48 | transition metal: semicore bands, d bands at E_F | rule | 12³ | 4³ | 8³ / 16³, 8³ |
| G4 | H3S, Im-3m, 200 GPa (F0c cell) | 4, 48 | high-pressure hydride: stiff modes, large bandwidth | 72 | 15³ | 3³ | 6³ / 16³, 8³ |
| G5 | LiAl2Rh (JVASP-18768, F0c cell) | 4, ≥ 24 | transition-metal compound with a soft lattice | 88 | 8³ | 2³ | 4³ / 10³, 5³ |
| G6 | TiN, rock salt, a 4.25 Å, **nr1 = nr2 = nr3 = 39 (3·13) forced by input** | 2, 48 | awkward FFT factor: outside QE's own allowed set and outside FDX | rule | 8³ | 2³ | 4³ / 12³, 6³ |
| G7 | MgB3C3__s1 (F0c re-relaxed cell) | 7, **2** | low symmetry: bases, SYMROT and STAR have little to share | 84 | 4 4 3 | 2 2 1 | 4 4 2 / 8³, 4³ |
| G8 | B13C2, R-3m, built as (B12)(CBC), vc-relaxed once before the runs | **15**, 12 | large cell: N_pw ≈ 10⁴, 45 modes, memory | rule | 4³ | 2³ | 4³ / 6³, 3³ |

Why each weak spot: the levers were tuned on 3- and 9-atom cells with 24 and 48 operations, N_pw 1566–5900, on two FFT grids, with no heavy d-band metal at E_F except Ir. G1–G8 each break one of these.
G2 is the optional 1-atom metal; it stays because it is cheap and it is the sharpest test of fixed per-launch costs.
G3 and G5 are the transition-metal pair; G5's cell comes from F0c as it stands.

**Fallback checks (not materials; ph.x only, `max` against all-off, frequencies compared, q2 on a small grid):**
- D1: bcc Fe, nspin = 2: magnetic, which the engine declines by design. Must run to the end with the vanilla loop and give the all-off answer.
- D2: MgO with `occupations = 'fixed'`: the true-insulator path of the same engine check.
- D3: MgB3C3__s1 at its funnel grade (k 6 6 7, q 2³): `max` to the "rtsc cbsr auto" line and through q2 only. This is the case that FUTURE-SPEED warned about (expected ≈ 10² bases of ≈ 0.6 GB), so it is the
  memory path on a real low-symmetry cell. It must print CG and must not exhaust memory.

Input preparation allowed before the runs, none of it comparing levers: B13C2's vc-relax, and one scf + nscf per material to read E_F and fix the frozen window, so all four sets share one epw.in.
If an input turns out to be unusable in preparation (for example a cell that will not converge), it is fixed there and the change is logged here before any lever run.

### Lever sets (switch names checked against the source at 17883ebc0; a flag is on only for the value `1`)
**all-off:** RTSC_STH=1 RTSC_STRICT=0 RTSC_STH_DVB=0 RTSC_W90_NT=1. (DVB and W90_NT are run-fork.sh defaults, so they are overridden, as in the scored all-off runs.)

**safe:** RTSC_STH=1 RTSC_STRICT=0 RTSC_STH_FFTDX=1 RTSC_STH_FDXDP=1 RTSC_STH_OBATCH=1 RTSC_STH_QCACHE=1 RTSC_NSCF_SYMROT=1 RTSC_PHQ_FUSE=1 RTSC_PHQ_VKB=1 RTSC_STH_DVB=1 RTSC_DVXC_GPU=1
RTSC_STH_KDEV=2 RTSC_STH_PSIR_PART=1 RTSC_STH_DRHO2=1 RTSC_PW_SYMROT=1 RTSC_W90_NT=8 RTSC_W90_WNT=8 RTSC_W90_MMN=1 RTSC_W90_AMN=1 RTSC_EPW_FDX=1 RTSC_EPW_BEC=1 RTSC_EPW_STAR=1 RTSC_EPW_B2W=1
RTSC_EPW_IDLEWAIT=1 RTSC_EPW_NOVME=1 RTSC_W90_NOHR=1.

**fast:** safe + RTSC_THRESH_L=10 RTSC_STH_DEFL=1.

**max:** safe + RTSC_STH_CBSR=1 RTSC_STH_CBSR_XCH=1 RTSC_STH_CBSR_HOST=1 with the automatic fallback on (the default), no RTSC_THRESH_L. XCH and HOST (class B) are included because the fallback's memory path is written for them.

Not in any set: EPW_NOWRITE (breaks EPW-only reuse), RUN_DIRECT, RUN_PREWARM, MP, DRF, TR2_SCALE.

### How a chain runs
run-fork.sh at the `gen` head, one binary (bin-gen) for every run, pw/ph 1 rank × 8 OMP, EPW 8 ranks + MPS, one chain at a time on the GPU, `flock /root/rtsc/TIMING.lock`, stages scf, ph, collect, nscf, epw.
Peak GPU memory: `nvidia-smi` sampled every 2 s during the chain (includes contexts and MPS). Each material runs all-off first, then safe, fast, max. Order: G1, G2, G3, G4, G6, G5, G7, G8, then D1–D3.
Spread: 3 repeats of G2's all-off ph.x, reported as (max − min)/mean. Speedups inside 1 + 3 × spread are not differences.
**A chain that runs over 2 × its all-off time + 300 s, or 45 min, is killed and recorded as a timeout** (all-off itself: 45 min). The cost line: if the projected total passes $5.40, the remaining work is cut in the
reverse of the order above, and the cut is reported.

Decision the long runs inform (RULES: runs over about 30 min): G7's and G8's all-off chains (≈ 30 min each, estimated) decide whether the safe, fast and max sets can be defaults on low-symmetry and large cells.

### Gate per material and set (the amended MILESTONE-1 gate, against that material's own all-off run)
- phonon frequencies: every mode at every q within max(0.5 %, 1 cm⁻¹); the three acoustic modes at Γ excluded;
- λ and ω_log within 1 %;
- Allen–Dynes T_c within max(2 %, 0.5 K);
- a true insulator (G1, D2): frequencies only.

λ, ω_log, T_c^AD are computed from prefix.a2f_proj with the funnel's `a2f_metrics` (EPW's formulas, first smearing column, μ* 0.10), so a zero-T_c case still gives numbers.
Margin used = the largest ratio of |deviation| to tolerance over all gated quantities, with the binding quantity named.

### Report per material and set
Speedup per stage and end to end; peak GPU memory; which levers declined and whether cleanly (a decline line, no crash, right answer); margin used.

### Verdict rule (fixed here)
- A lever set is fit to be a DEFAULT only if **every** material passes the gate with no crash.
- A failure restricts that lever to the domain where it passed. The failure stays in the table.
- A crash is any stage exiting non-zero, or a timeout under the rule above. A clean decline is not a failure.
- Report the median and the worst speedup, plainly, across all materials that ran.
- If a material's all-off run itself fails, its row stays in the table marked "no reference". That material cannot pass, and the verdict says so.
- Wannier: spread-sum outliers and stuck runs are noted, not fixed (another agent's work). A gate failure that traces to Wannier stays a failure in the table; a diagnostic may sit next to it, labelled as such, and does
  not change the verdict.
- A set failing because of one lever may be followed by a diagnostic to name that lever. The diagnostic is labelled, and it changes the verdict only in naming the domain, as the second line above says.

### ETA (written 02:30 UTC)
Part 1 finish ≈ 03:15. Preparation ≈ 03:45. Chains (estimates, all-off / each lever set): G1 100 / 60 s, G2 200 / 120, G3 300 / 180, G4 400 / 250, G6 300 / 180, G5 1200 / 700, G7 1800 / 900, G8 2000 / 1100; D1–D3 ≈ 25 min.
Sum ≈ 5 h, so the end is ≈ 09:00 UTC. Total ≈ 7 h ≈ $4.1, plus the save. Cap $6 (10 h).
