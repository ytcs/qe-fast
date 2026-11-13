# Changelog

All notable changes to qe-fast are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

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

### Removed

- Switches that were measured and not adopted: `QEF_THRESH_G` (now fixed at 100), `QEF_MIX_Q0` (Kerker mixing,
  `LR_Modules/qef_mix.f90`), `QEF_NSCF_ETHR_SCALE`, `QEF_W90_NITER` and `QEF_MP_INPLACE`. `UtilXlib/mp_base.f90` is
  upstream again.
- The `./configure` and `make` build, dev-tools, archive and the upstream CI files. CMake is the only build.
- The upstream packages CPV, COUPLE, GUI, GWW, KCW, PWCOND, QEHeat and XSpectra, their test-suite cases, and the
  pw2qmcpack plugin.

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
