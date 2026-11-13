# qefast/

Everything specific to qe-fast that is not a change to an upstream file lives here. See the
[top-level README](../README.md) for what the project is and how to build it.

- `SWITCHES.md`: every run-time switch, with its default, class, gate evidence and fallback.
- `config/`: `qef_config`, the central preset table and `qef_getenv`.
- `sth/`: the batched Sternheimer and nscf engines.
- `sth/fdx/`: the cuFFTDx kernels, a stub for builds without MathDx, and `deps.sh`.
- `epw/`: GPU paths for EPW, and the gate and comparison scripts.
- `w90/`: patched copies of two wannier90 routines and the scripts that generate them.
- `pod/`: setup and build scripts for an NVHPC container.
- `pseudo.sh`: fetches the SG15 pseudopotentials that the benchmarks use.
- `examples/eph-tc/bench/`: benchmark inputs (`B1` is MgB₂, `B2` is Mg₂IrH₆, each `smoke` and `frozen`), the chain
  runner `run-chain.sh`, and `gen/`, the held-out study's inputs and gate report.
- `examples/eph-tc/funnel/`: a tiered candidate-screening pipeline built on these tools, kept as an example.
  It needs Python packages and data that are not part of this repository (see NOTICE).

## Branches and tags

- `main`: integration branch.
- `qe-7.6`: the upstream release tag this fork branched from.
- `v*`: qe-fast releases.
