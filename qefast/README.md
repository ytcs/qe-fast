# qefast/

Everything specific to qe-fast that is not a change to an upstream file lives here. See the
[top-level README](../README.md) for what the project is and how to build it.

- `SWITCHES.md`: every run-time switch, with its default, class, gate evidence and fallback.
- `config/`: `qef_config`, the central preset table and `qef_getenv`.
- `sth/`: the batched Sternheimer and nscf engines.
- `sth/fdx/`: the cuFFTDx kernels, a stub for builds without MathDx, and `deps.sh`.
- `epw/`: GPU paths for EPW, and the gate and comparison scripts.
- `w90/`: patched copies of two wannier90 routines and the scripts that generate them.
- `gw/`: the GW stages as the library `qe_gw` (pair densities, χ₀ by Kramers–Kronig, ε⁻¹, Im W; host and GPU),
  `chi0.x`, and `chi0_ref.py`, the Python reference they are checked against.
- `lint/`: source checks run in CI.
- `compat/`: compatibility shims for the plugins.
- `pod/`: setup and build scripts for an NVHPC container.
- `pseudo.sh`: fetches the SG15 pseudopotentials that the benchmarks use.
- `examples/eph-tc/bench/`: benchmark inputs (`B1` is MgB₂, `B2` is Mg₂IrH₆, each `smoke` and `frozen`), the chain
  runner `run-chain.sh`, and `gen/`, the held-out study's inputs and gate report.

## Scope

qe-fast is a general ab initio package. A result belongs here when it depends only on the crystal and the level of
theory and can be checked against another code or an exact identity: bands, phonons, electron–phonon matrix
elements, ε⁻¹_GG′(q, ω). Anything that adds a probe, a scenario or a ranking, or whose acceptance test is a number
from an application, lives in an application repository. Such repositories pin a qe-fast commit and build their
programs against its libraries with `-DQEF_EXTRA_DIRS=<dir>` (top-level `CMakeLists.txt`).

## Branches and tags

- `main`: integration branch.
- `qe-7.6`: the upstream release tag this fork branched from.
- `v*`: qe-fast releases.
