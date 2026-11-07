# Security policy

## Reporting a vulnerability

Report vulnerabilities privately through GitHub: on the repository's **Security** tab, choose
**Report a vulnerability**. Please do not open a public issue for them.

Reports are handled on a best-effort basis. This is a small volunteer project with no response-time commitment.

## Supported versions

Only the latest release and `main` get fixes.

## Scope

In scope: code that qe-fast added or changed (`qefast/` and the files listed in `NOTICE`) and its scripts.

A problem in unmodified Quantum ESPRESSO, EPW, wannier90 or another dependency belongs to that project.

qe-fast is a scientific code. It trusts its input files and run environment and is not designed to process untrusted
input. The scripts under `qefast/pod/`, `qefast/pseudo.sh` and `qefast/sth/fdx/deps.sh` download third-party software
and data, and `qefast/pseudo.sh` does so over plain HTTP with no checksum, so run them only on hosts where that is
acceptable.
