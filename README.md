# qe-fast

qe-fast is a fork of [Quantum ESPRESSO](https://www.quantum-espresso.org) 7.6, EPW included. It adds a
batched GPU path for the phonon and electron-phonon stages of a superconductivity workflow: `pw.x` scf,
`ph.x`, `pw.x` nscf, wannier90 and `epw.x`.

It is for people who run DFPT and EPW on NVIDIA GPUs for non-magnetic, norm-conserving metals and who are
prepared to check what a run tells them. Outside that scope qe-fast runs Quantum ESPRESSO's own code, with no
speedup, so upstream QE is the better choice there.

This is version 0.1.1. The upstream README is kept as [README.QE.md](README.QE.md); the upstream documentation
in `Doc/` and `*/Doc/` applies to everything that qe-fast did not change.

Upstream packages that neither of qe-fast's uses needs have been removed: CPV, COUPLE, GUI, GWW, KCW, PWCOND,
QEHeat, XSpectra and the pw2qmcpack plugin, together with the `./configure` and `make` build. NOTICE lists them.

> **Not affiliated.** qe-fast is an independent project. It is not affiliated with or endorsed by the Quantum
> ESPRESSO Foundation or the EPW developers. The Quantum ESPRESSO name is used descriptively only, to say what
> this software was derived from. Please do not send qe-fast questions to the upstream mailing lists or issue
> trackers.

## Build

CMake only; the `./configure` and `make` build has been removed. Clone with git:
CMake fetches the submodules, so GitHub's source archives will not build.

    git clone https://github.com/ytcs/qe-fast && cd qe-fast

### GPU

Needs the NVIDIA HPC SDK. `qefast/pod/build-gpu.sh` is written for NVHPC 26.5 at its default install path in
a Linux container; edit the `NV` line for another location. It is the only toolchain that has been tested.
The argument is the CUDA compute capability (86 = A40, 90 = H100, the two GPUs the results below come from).
`CPU=0` skips the script's additional CPU build.

    CPU=0 qefast/pod/build-gpu.sh . 90          # binaries in build-gpu/bin

That build contains no MathDx code. The cuFFTDx kernels (`QEF_STH_FFTDX`) compile to a stub that always
declines, and the run keeps cuFFT. To compile the real kernels, for your own use only (see
[Licensing](#licensing)):

    export FDX_DEPS=$HOME/.cache/qefast-fdxdeps # the script's default is under /root
    eval "$(qefast/sth/fdx/deps.sh)"             # fetches nvidia-mathdx 25.6.0 and CUTLASS 3.9.2
    CPU=0 qefast/pod/build-gpu.sh . 90

Running `deps.sh` means accepting NVIDIA's licence for MathDx. The cmake output says which one you got: `qef fftdx: cuFFTDx kernels for sm_90` or `qef fftdx: stub`.

### CPU

    sudo apt-get install gfortran gcc cmake ninja-build libopenblas-dev libfftw3-dev
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DQE_ENABLE_MPI=OFF -DQE_FFTW_VENDOR=FFTW3
    cmake --build build --target pw.x ph.x epw.x

`pw.x`, `ph.x` and `epw.x` build on CPU, but the GPU engine isn't compiled in, so there is no speedup to expect.
CI builds them without MPI and runs five `test-suite` cases against the upstream benchmarks (pw, ph, epw, tddfpt,
hp). No preset has been gated on a CPU build. For CPU-only work, use upstream QE.

## Quick start

The programs take the same inputs as QE's. Run `pw.x` and `ph.x` as one MPI rank per GPU with OpenMP threads,
and pass any `QEF_*` variables on to remote ranks. Nothing needs to be set: each lever checks the input, the build
and the device and turns itself on only where it applies, then prints one line saying so, for example

    qe-fast: automatic; each lever reports below whether it applies to this run
    qe-fast: GPU Sternheimer engine declined: not a smeared metal (fixed occupations or tetrahedra)
    qe-fast: ph.x nscf by symmetry rotation on

The automatic set is the `safe` preset below: only levers that keep the output bitwise identical or differ at
rounding level. A lever that declines leaves QE's own code path in charge.

The MgB₂ smoke chain runs every stage, in about 430 s under `safe` on an A40. It needs `mpirun` on the `PATH`
(NVHPC ships one in `comm_libs/mpi/bin`):

    qefast/pseudo.sh                            # fetches the SG15 pseudopotentials into qefast/pseudo/sg15
    QEF_PRESET=safe qefast/examples/eph-tc/bench/run-chain.sh \
        qefast/examples/eph-tc/bench/B1/smoke build-gpu/bin /tmp/b1-safe 8 8

The arguments after the work directory are the OpenMP threads for `pw.x` and `ph.x` and the MPI ranks for
`epw.x`. Stage times go to `timings.tsv` in the work directory.

## Presets

Presets are a developer setting, for benchmarks and gates; leave `QEF_PRESET` unset for normal runs.
`QEF_PRESET` selects a set of switches. Any `QEF_*` variable you set overrides its preset entry. The full table
of switches, with the class, evidence and fallback of each, is in [qefast/SWITCHES.md](qefast/SWITCHES.md).

| preset | what it is | status |
|---|---|---|
| `off` | The batched engine with every optional lever off. This is not plain QE: `QEF_STH=0` gives QE's own loop. | Reference for every gate |
| `safe` (= automatic, the default) | Levers that keep the output bitwise identical or differ at rounding level, and leave the Wannier gauge alone | Full gate passed on 4 held-out materials and on the two smoke benchmarks |
| `fast` | `safe` plus lever C and three Wannier-related switches | **Phonon-gated only** on held-out materials |
| `max` | `safe` plus the CBSR exact solve and the same three switches | **Phonon-gated only** on held-out materials |

The gate, against the same material's `off` run: every mode at every q within max(0.5 %, 1 cm⁻¹), the three
acoustic modes at Γ excluded; λ and ω_log within 1 %; Allen–Dynes T_c within max(2 %, 0.5 K).

`fast` and `max` passed the frequency gate on all eight held-out materials, but they may land on a different
Wannier solution from `off`. On four of the eight, λ, ω_log or T_c failed the gate (λ moved by up to 29 %, on
TiN), and B₁₃C₂ has no λ reference. `off` is not a truth reference for λ, but the result is not the same
answer. Do not use `fast` or `max` for λ or T_c on a new material without checking it against `off`.

## Supported cases and fallbacks

The batched `ph.x` engine runs when all of these hold; `LR_Modules/sth_batched.f90` is the authority:

- norm-conserving pseudopotentials;
- collinear and non-magnetic (`nspin_mag = 1`), no spin-orbit coupling;
- smearing occupations, so no `occupations = 'fixed'`, no tetrahedra and no two-chemical-potential runs;
- no DFT+U, meta-GGA, EXX, `real_space`, scissor or `lelfield`;
- k-point sampling other than Γ-only (`gamma_only`);
- serial FFT within a pool, one band group, no task groups. Any number of pools.

Any other case takes QE's own loop automatically, with no speedup. That includes magnetic, SOC, USPP/PAW and
DFT+U systems. `QEF_STRICT=1` turns such a fallback into a stop.

`QEF_STH_DVB=1` falls back to QE's own dV_bare loop when the engine declines at run time, and with
`QEF_EPW_NOVME=1` a run-time decline of the fused fine path hands the run to the stock fine path. Before 2026-11-13
both could stop a run; the new fallbacks have not yet run on a GPU.

## Performance

Measured numbers only. **Every figure below is relative to qe-fast's own code, not to unmodified Quantum
ESPRESSO.** There is no end-to-end comparison against unmodified QE on the same hardware, because the
full-size unmodified runs were not completed, so qe-fast makes no such claim.

### Tuned benchmarks, H100

The two benchmarks that the code was tuned on, at their full size ("frozen" inputs), scf through EPW. Each run
used a hand-picked lever set, not a preset. Neither result is what the default gives.

| system | baseline | this run | speedup | gate against the baseline |
|---|---|---|---|---|
| MgB₂ | 2849 s | 240.84 s | **11.83×** | passed |
| Mg₂IrH₆ | 4030 s | 818.08 s | **4.93×** | passed; **the 10× target was not reached** |

- Hardware: one H100 80GB HBM3. `pw.x` and `ph.x` as 1 rank × 8 OpenMP threads, `epw.x` as 12 ranks with MPS.
- Baseline: qe-fast's first complete engine (commit `94c1e2ccb`, 2026-10-28) on an H100, end to end. For MgB₂,
  an all-levers-off run of the same code on the same host took 2562.21 s, which gives 10.64×.
- Both runs were at commit `faba36128` (2026-10-30), and they have not been repeated on the v0.1.0 tree.
- Both used `QEF_STRICT=1`, the cuFFTDx kernels, and `QEF_W90_WNT` and `QEF_W90_MMN`, which `safe` leaves out. A
  build without MathDx cannot reproduce them. The kernels are compiled for the FFT grids of these two systems.
- MgB₂ used the CBSR solve, `QEF_STH_DEFL=1` and `QEF_EPW_NOWRITE=1`, without `QEF_PW_SYMROT`.
- Mg₂IrH₆ ran without CBSR, with `QEF_THRESH_L=10` and `QEF_PW_SYMROT=1`.

### Default preset on held-out materials, A40

End to end, `safe` against `off`: same A40, same inputs, same material. `off` is the engine with its levers off,
not plain QE. The four held-out materials are the ones that failed an earlier, larger `safe`, so they are the
hardest tests of the accuracy gate and not a random sample.

| | Al | Nb | LiAl₂Rh | TiN (39³ grid forced) |
|---|---|---|---|---|
| `safe` end to end | 1.50× | 1.35× | 1.07× | 1.10× |

The median is 1.23×, the range 1.07–1.50×. All four passed the full gate. `safe` itself was not run on the other
four held-out materials (MgO, H₃S at 200 GPa, MgB₃C₃, B₁₃C₂). On the two tuning systems at smoke size, on the
same A40, `safe` gave 1.50× (MgB₂) and 1.61× (Mg₂IrH₆), also passing the full gate.

### `fast` and `max` on held-out materials, A40

End to end against `off`, median / worst over the seven held-out materials with a complete chain: `fast`
1.68× / 1.42×, `max` 1.95× / 1.09×. These two presets are phonon-gated only, see [Presets](#presets).

## Limitations

- **FFTDX is compiled for two FFT grids only**: 30×30×36 and 54×54×54. It declined on all eight held-out
  materials, which is much of why `safe` gains little there: `safe` left the solve stage where it was on five of
  seven. The decline is clean and the run keeps cuFFT.
- **No speedup for magnetic, SOC, USPP/PAW or DFT+U systems**, or for the other cases listed above. They take
  QE's own loop.
- **`fast` and `max` are phonon-gated only** on held-out materials and may land on a different Wannier solution.
- **λ is sensitive to the Wannier gauge.** Rounding-level changes can move a run off a saddle in the spread, and
  moving Wannier functions by lattice vectors changed λ by 3.7 % on a Mg₂IrH₆ test at unchanged spread. When you
  compare runs, compare Cartesian Wannier centres, not spreads
  (`qefast/examples/eph-tc/bench/gen/wfcentres.py`).
- **Evidence is limited**: two tuning systems and eight held-out materials, on A40 and H100 GPUs, with SG15
  pseudopotentials. Gate runs were GPU runs. The grids of the held-out runs were trimmed for cost.
- One toolchain has been tested for GPU builds (NVHPC 26.5).
- CPU builds of `ph.x` have been run on one small test (Si at X, serial and MPI); `epw.x` on CPU has not been run.
- CI covers a CPU build and a small `pw.x` test. It does not run GPU code, so the GPU results above are not
  regression-tested on each commit.

## Licensing

qe-fast is free software under the GNU General Public License, version 2 or (at your option) any later version,
as is Quantum ESPRESSO. The licence text is in [License](License); [NOTICE](NOTICE) records the modifications
and every third-party component.

**Do not distribute binaries built with MathDx/cuFFTDx.** NVIDIA's licence for it is proprietary and forbids use
that would subject it to an open-source licence, while the GPL would require exactly that of a distributed
binary. Nothing of MathDx is in this repository, and builds for your own use are not affected. For a binary you
pass on, leave `QEF_MATHDX_INC` unset. See NOTICE for the full position. This project publishes no binaries.

## Citing

If you use qe-fast in published work, cite it and the upstream Quantum ESPRESSO and EPW papers. They are listed
in [CITATION.cff](CITATION.cff).

## Project files

[CHANGELOG.md](CHANGELOG.md), [CONTRIBUTING.md](CONTRIBUTING.md), [SECURITY.md](SECURITY.md),
[qefast/README.md](qefast/README.md) (layout of the fork).
