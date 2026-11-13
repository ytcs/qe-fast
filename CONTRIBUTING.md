# Contributing to qe-fast

qe-fast is an independent fork of Quantum ESPRESSO. Use this repository's
[issues](https://github.com/ytcs/qe-fast/issues) and pull requests for anything specific to it. Do not send qe-fast
reports to the upstream mailing lists or trackers. Report a problem upstream only if you can reproduce it with unmodified
Quantum ESPRESSO, through [their channels](https://www.quantum-espresso.org). Report security problems as described in
[SECURITY.md](SECURITY.md).

## Bug reports

Give the commit, the GPU and toolchain, the `qe-fast` lines at the top of the output and any `QEF_*` variables, and a small input that
shows the problem. If the problem is a wrong number, say whether it also happens with `QEF_PRESET=off` and with
`QEF_STH=0`.

## Pull requests

- Open them against `main`. CI must pass: a CPU build, a small `pw.x` test, and syntax and naming checks.
- A change that can touch results must state its class (see [qefast/SWITCHES.md](qefast/SWITCHES.md): B bitwise, E
  rounding-level, N changes the numbers) and give the gate evidence: GPU runs against the same material's `off` run.
  CI does not run GPU code, so that evidence has to come from you.
- New run-time switches go through `qef_getenv` and get a row in `SWITCHES.md`. A new lever is off in every preset until
  it has passed the gate.
- Keep speed claims measured, with hardware, baseline and scope. Say what did not improve.
- Prefer new files under `qefast/`, and keep changes to upstream files small and marked with a `qe-fast:` comment.
- Comments describe the code. They do not point to notes that are not in this repository.

## Licensing of contributions

By contributing you agree that your contribution is distributed under GPL-2.0-or-later. New source files start with:

    Copyright (C) <year> The qe-fast authors
    SPDX-License-Identifier: GPL-2.0-or-later

When you change an upstream file for the first time, add it to the table in `NOTICE` and add a `qe-fast: modified` line
at its top. Keep every existing copyright notice.

Before importing third-party code, check its licence and record it in `NOTICE`:

| Package | Licence | Status |
|---|---|---|
| wannier90 (pinned commit) | GPL v2 | compatible |
| ABINIT | GPL v3 | do not merge code; a GPL-v3-only part would remove the v2 option; use it as an external cross-check |
| python-sscha / CellConstructor | check before import | not imported |
| phonopy, pymatgen, ASE | BSD / MIT / LGPL | interface-level use only |

Never commit MathDx or cuFFTDx files, or anything else under a proprietary licence. See `NOTICE`, section 3.
