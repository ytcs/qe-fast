# Changelog

All notable changes to qe-fast are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `chi0.x` (`qefast/gw/chi0.f90`, issue #1): the inverse dielectric matrix ε⁻¹_GG′(q, ω) by sum over bands with
  local fields (RPA), for norm-conserving insulators; Kramers–Kronig by FFT on the energy grid; full-matrix output.
  CUDA builds run it on the GPU (`chi0_dev.f90`): Si 8³, 60 bands, 3000 energies in 5.4 s on an H100. Both paths equal
  the Python reference `qefast/gw/chi0_ref.py` element-wise to 2×10⁻¹³; the reference agrees with `turbo_eels` on Si
  (static ε within 0.3 %, loss peaks within 0.1 eV). CI checks the host path against the reference.
- `qefast/gw/gamma_ref.py`: reference contraction of Im W with pair densities for the impact-ionization rate (E-13).
- `gamma_ii.x` (`qefast/gw/gamma_ii.f90`, issue #1): the impact-ionisation rate Γ_ii = 2|Im Σ|/ħ of every hot
  electron up to a carrier cutoff, from the RPA Im W of every q on the grid, and the partial rate into every final
  state (for the cascade). Host path; on Si 4³ it equals `gamma_ref.py` to the printed precision (7 digits) and takes
  half its time. The q = 0 term is left out (weight 1/N_k). It shares its stages with `chi0.x` through
  `qefast/gw/gw_common.f90`. CI checks it against the reference.

### Fixed

- Ultrasoft Q(q) for q ≠ 0 (`upflib/qvan2.f90`, `compute_qqr`/`compute_qqc`): |q|² was q₁² + q₂² + q₁², so the
  third component never entered. This affected EELS (`turbo_eels.x`) and Berry-phase/finite-field runs with ultrasoft
  pseudopotentials on every build. Q(q) was also left unset for norm-conserving species. Backported from upstream
  development (QE commits b987d733 and 24c9753d, listed for 7.6.1). On `tddfpt_eels-al` (Al_uspp) the CPU ε⁻¹ moves
  by up to 2×10⁻³; the case now equals upstream's regenerated benchmark (taken from QE commit 4aead329) to 3×10⁻⁵.
  `pw_berry` and `pw_electric` still pass.
- `turbo_eels.x` with ultrasoft pseudopotentials on GPU, Sternheimer solver: ε⁻¹ was off by 2–7 % (Al_uspp on an
  H100, the same in vanilla 7.6). `orthogonalize_omega` built `becp` on the host while `s_psi` read the stale device
  copy; `ccg_psi`, the preconditioner of the frequency-dependent solver, ran on host copies the solver never used.
  Ultrasoft with and without spin–orbit, collinear and noncollinear, now equal the CPU on the H100 to 4×10⁻⁶ (ε⁻¹ at
  1.44 eV); the norm-conserving case was right before and is unchanged. The same upstream commit (QE 4aead329) fixes
  host/device mixing in `lr_sm1_psi`, used by the Lanczos solver.
- Known issue: ultrasoft Lanczos EELS on GPU still differs from the CPU (Al, collinear: β₁ 0.5962 against 0.5940
  after the fixes below; norm-conserving agrees to 10⁻⁶). Being traced.
- TDDFPT on GPU with ultrasoft pseudopotentials: several routines filled `becp` with the host `calbec` and then
  called the offloaded `s_psi`, which read the stale device copy. Affected: `lr_apply_s`, which applies S to every new
  Lanczos vector, plus `sd0psi`, `lr_read_wf`, `lr_dvpsi_e` and `lr_normalise`. `becp` now goes to the device in
  between. `lr_sm1_psi` built the k+q projectors on the device only and then used the host copy. Not yet measured.
- EPW iterative BTE in non-MPI builds: `transport_ibte` opened the merged sparse-matrix files without their image
  suffix, created empty files and stopped ('error in reading sparse_q'). Upstream development has the same lines.
- XClib, LAXlib and upflib error routines flush the message before `MPI_Abort`. On MPI builds an input error such as
  SCAN without Libxc used to end with an empty output.
- Every bug fix upstream QE has made since the 7.6 release, about 90 commits from the development branch, is
  backported. Each is cherry-picked with its upstream author and commit id; the commit history lists them. Among
  those that change results:
  - hybrid-functional stress (the q = 0 term was missing, the normalisation was wrong) and the Gau-PBE density
    derivative;
  - averaging of fully relativistic pseudopotentials used without spin–orbit;
  - meta-GGA: the stress for ecutrho > 4 ecutwfc, Libxc's sigma floor, and the density thresholds;
  - magnetic noncollinear phonon symmetrisation;
  - tetrahedra with fixed magnetisation under MPI, and eigenvectors of the "direct" LAXlib diagonalisation;
  - BFGS (s·y after Powell damping, trust radius in the line search), velocity-Verlet MD, and the vc-relax charge
    extrapolation;
  - DFT+U: noncollinear symmetrisation, orbital-resolved U, `read_ns`, and HP with one k-point per pool;
  - epsilon.x (χ was 4× too large for nspin = 1), the pw2wannier90 `CONJG` and spinor factors, and pw2wannier90's
    hydrogenic radial integral (up to 99 % off for r = 3);
  - `adddvscf` called with two of its three arguments (EELS and el-ph with ultrasoft pseudopotentials).

  GPU data-synchronisation fixes cover CG eigenvalues, `vloc_psi`, the Hubbard forces, PAW noncollinear and
  `doverlap_inv`. Crash and hang fixes cover species labels longer than 3 characters, `read_rhog`, `wfcinit`,
  concurrent runs in one directory, and re-initialisation after `clean_pw`. Two local corrections were needed:
  - `rotate_and_add_dyn` keeps its 7.6 interface for the d3q plugin, which calls it without the new time-reversal
    flags;
  - `q2qstar.x` passes `t_rev` in the right position (upstream passes it after `rtau`).

  Benchmarks were taken from upstream only where our build reproduces them: the hybrid, meta-GGA, relax, constraint,
  tddfpt_hyb, pw_irrbz and EPW sets. `pw_workflow_vc-relax_dos` is regenerated here, since the BFGS fix changes the
  path and both endpoints are converged. On the serial CPU build, 74 of the 75 pw, ph, hp, pp and tddfpt test
  directories pass against the updated benchmarks; `ph_multipole` needs `spglib` in the environment. 19 of the 21
  EPW directories pass; `epw_lsda_metal` and `epw_mob_U` need pseudopotentials the QE server does not host.
  `tddfpt_hyb` failed on vanilla 7.6 too and now passes. On an H100 the 32 GPU-relevant test directories pass,
  including `ph_metal`, which vanilla 7.6 fails there; the one exception is `pw_metaGGA`'s SCAN case on a build
  without Libxc.
- ph.x, `electron_phonon='interpolated'`: the linewidths depended on the basis the nscf diagonalisation returned
  within degenerate sets, and on the eigenvector basis of degenerate phonon modes. Upstream has the same dependence.
  On `ph_metal` (Al, q at X) γ moved by up to 8 % with the nscf starting vectors, even with both starts converged to
  ethr = 3×10⁻¹⁴. |g|² is now averaged over the degenerate bands at k and k+q before interpolation, and γ, λ over
  frequency-degenerate modes. The two starts agree to the printed 0.01 GHz and the benchmark to 0.57 GHz.
- EPW `pw2wan`: hydrogenic radial integrals. The logarithmic mesh cut at 10/α truncated the r = 2 and r = 3 guiding
  functions, the bug upstream fixed in pw2wannier90 but not in EPW's copy. The r = 3 function also had the wrong
  normalisation exponent (2/3 for 3/2) and linear coefficient (1.5 for 2/3). Both radial routines are fixed.

- XML read of `pseudo_file`: the `l` attribute was read and broadcast without being allocated. FoX builds crashed on
  every read of a saved run (`ph.x`, `gipaw.x`, ...); other builds silently dropped it.
- The qe-gipaw plugin builds and runs: pinned at upstream 02e76b6 (7.6 tag plus a format fix). Its quartz NMR test
  reproduces the upstream 7.6 reference shieldings.
- The d3q plugin links against QE 7.6 through `qefast/compat/d3q_compat.f90`.
- FoX configures with CMake 4.
- `ph.x` stopped in every CPU build ("howmany different from 1, not yet implemented for FFTW3", or "not yet
  implemented for parallel driver" with more than one rank per pool): the batched FFT in `dvqpsi_us` was taken
  without checking the build. It is now used in CUDA builds only. Si at X (`test-suite/ph_base`) runs serial, with
  pools and with a distributed FFT.
- `epw.x` stopped in every CPU `ephwrite` + Eliashberg run under the default `safe` preset
  ("QEF_EPW_NOVME=1 needs the fused fine path"): the switch assumed the GPU fine path. It now turns itself off in
  CPU builds and when the input rules out the fused path. Found by `test-suite/epw_super`.
- `QEF_PW_SYMROT` (`fast`) stopped a `pw.x` nscf run with `disk_io = 'none'`; it now declines and every k-point is
  diagonalized. `QEF_NSCF_SYMROT` (`safe`, the nscf inside `ph.x`) shares the same decline function now; a stop
  there was not observed.
- GPU (nvfortran) builds stopped at `wannier90.x`: `qef_config.mod` was not visible to the users of the patched
  wannier90 modules.
- Upstream QE 7.6 bugs in serial (`-DQE_ENABLE_MPI=OFF`) EPW builds, found by CI: `epwread` read the
  Wannier e-ph vertex into an unallocated `epmatwp` (segfault), and `etf_mem = 1` opened `epmatwp` with records
  longer than the slices it reads and read past the end of the file. `test-suite/epw_base` failed 7 of 15 serially on
  vanilla 7.6 and on the fork; it passes 15 of 15 now.
- An `epwread` run on files written under `QEF_EPW_NOVME=1` (a GPU `safe` Eliashberg run) crashed reading the
  empty `vmedata.fmt` (`test-suite/epw_super`, GPU). The files now carry a marker line; the reader leaves the
  velocity matrices at zero and stops only if the new run needs them.
- The `dvxc_gpu` startup line no longer says "on the device" in builds without OpenACC.

### Changed

- The GPU Sternheimer engine runs insulators with fixed occupations (the projector onto the occupied bands of k+q,
  as QE's `orthogonalize`). H100, against vanilla 7.6 on the same GPU: Si `ph.x` 501 s → 12 s, MgO 78 s → 6 s,
  frequencies within 0.02 cm⁻¹, ε∞ and Born charges equal in every printed digit (the electric-field solve ran on the
  engine); the GPU `ph` test category has the same 87 passes and the same two failures as vanilla: ph_raman with the
  same values, ph_metal's λ entries 34.06 / 34.95 against vanilla's 34.10 / 35.00 (tolerance 1.5 from the benchmark).
- Automatic configuration: with `QEF_PRESET` unset (the default, `safe`), the run's output says so and the main
  levers print one line each, on or declined with the reason (`qef_note`). The GPU engine's decline (an
  insulator, USPP/PAW, spin-orbit, ...) is now visible in the output.
- `QEF_STH_DVB` falls back to the stock dV_bare loop when the engine declines at run time (device memory) instead of
  stopping the run.
- `QEF_EPW_NOVME`: when the fused fine path declines at run time (pair-kernel shared memory, images), the stock
  fine path now runs instead of stopping; neither path reads the skipped velocity matrices without a velocity
  consumer. It declines up front with quadrupoles, which need the Berry connection.
- `QEF_STH_KDEV`'s device cache must leave 4 GiB free, else it is dropped; vkb chunks are capped at 512 MiB.
  `QEF_STH_PSIR_PART` reserves room for the NL state and FDXDP's index map allocated after it.
- `ph.x` recover with bands read from disk forgets the SYMROT plan of an earlier q.
- After a dV_bare fallback the engine can come back later in the same q: an irrep that asks for dV_bare on the
  engine rebuilds the per-q state when it was built without it (an illegal address on the H100 forced-decline test).
- The batched nscf engine (`QEF_STH=1`) converges max(2, nbnd/10) buffer bands above the ones it returns. Without them
  it returned the highest band of `test-suite/ph_ahc_diam` 1.15 eV too high at one k (36.06 against 34.91 eV; the
  lower Ritz value is the better one), and AHC's upper Fan was off by up to 2.3 % on GPU under every preset. With
  them ph_ahc_diam passes 16 of 16 on an H100, and B1 smoke stays inside the gate against `off` (0.069 cm⁻¹; 39 s).
  The run reports the largest residual of the returned bands. (For a while AHC runs declined the batched nscf.)
- `phq_fuse_all` sizes its k-point chunks from the free device memory (at most 1.5 GB, half of what is free above
  512 MiB); with 2 GiB free it ran out of memory at the second q-point.
- `QEF_W90_NT` never exceeds `omp_get_num_procs()`, and stays serial with `timing_level > 1`.

### Removed

- Switches that were measured and not adopted: `QEF_THRESH_G` (now fixed at 100), `QEF_MIX_Q0` (Kerker mixing,
  `LR_Modules/qef_mix.f90`), `QEF_NSCF_ETHR_SCALE`, `QEF_W90_NITER` and `QEF_MP_INPLACE`. `UtilXlib/mp_base.f90` is
  upstream again.
- The `./configure` and `make` build, dev-tools, archive and the upstream CI files. CMake is the only build.
- The upstream packages CPV, COUPLE, GUI, GWW, KCW, PWCOND, QEHeat and XSpectra, their test-suite cases, and the
  pw2qmcpack plugin.
- Development scaffolding: the offline engine harness and its dump hook (`QEF_STH_DUMP`), the FFT and cuFFTDx
  microbenchmarks, and the pod cache scripts.

## [0.1.1] - 2026-11-06

### Fixed

- `ph.x` and `epw.x` now compile without CUDA. GPU builds are unchanged. CI builds all three executables.

## [0.1.0] - 2026-11-06

First public release, based on Quantum ESPRESSO 7.6 (tag `qe-7.6`), EPW included.

### Added

- Batched GPU Sternheimer engine for `ph.x`, with an automatic fallback to QE's own loop for unsupported cases.
- Batched nscf, with symmetry rotation of wavefunctions in place of per-q eigensolves.
- GPU paths for the EPW coarse, fine and Eliashberg stages; OpenMP over k in two wannier90 routines.
- Optional cuFFTDx local operator for two FFT grids (30×30×36 and 54×54×54), with a stub for builds without MathDx.
- Optional CBSR exact solve with an automatic memory and payoff check.
- Central switch table `qef_config` with the presets `off`, `safe`, `fast` and `max`, documented in
  `qefast/SWITCHES.md`. The built-in default is `safe`.
- Benchmark inputs, chain runner, accuracy-gate scripts and held-out study inputs under `qefast/examples/eph-tc/`.
- README, NOTICE, CHANGELOG, CITATION.cff, SECURITY.md, issue and pull request templates, and CPU CI.

### Changed

- `external/devxlib` is vendored at `a6b89ef` instead of being a submodule, and carries the upstream MIT licence text.
- The upstream README moved to `README.QE.md`; `CONTRIBUTING.md` was rewritten for this fork.

### Known limitations

- Only CMake builds are supported, and only NVHPC 26.5 has been tested for GPU builds.
- CPU builds provide `pw.x` only: `ph.x` and `epw.x` do not compile without CUDA.
- `fast` and `max` are phonon-gated only on held-out materials and may land on a different Wannier solution.
- Magnetic, SOC, USPP/PAW and DFT+U systems take QE's own loop, with no speedup.
- `safe` gave 1.07–1.50× end to end against `off` on four held-out materials (A40).
- The H100 results on MgB₂ (11.83×) and Mg₂IrH₆ (4.93×, below the 10× target) used hand-picked lever sets, are relative
  to an earlier qe-fast build and not to unmodified QE, and have not been repeated on this tree.
- Binaries built with MathDx must not be distributed. See `NOTICE`.

[0.1.1]: https://github.com/ytcs/qe-fast/releases/tag/v0.1.1
[0.1.0]: https://github.com/ytcs/qe-fast/releases/tag/v0.1.0
