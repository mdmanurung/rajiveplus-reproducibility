# RaJIVE+ paper workflow: stage-0 foundation

**Pilot infrastructure only. Scientific readiness: NOT_READY. All methods,
simulations, calibration, clinical inference and scientific figures: NOT_RUN.**
This independently authored repository does not vendor or modify algorithms.
The public paper remote was checked with `git ls-remote` on 2026-10-22 and
returned no refs. Development is local only; no pushes or remote changes.

## Executable scope

Seven exact-source processed COVID CSVs → byte-checksummed acquisition receipt
→ defensive audit → three local sample manifests → source/environment/run
receipt. Snakemake runs locally, one core. SLURM is unavailable and unused.
Input size is 44,613,388 bytes (seven files); no method fitting occurs.
See `config/datasets.yaml` for the full source SHA, paths and seven SHA256s.

```bash
uv venv .venv --python 3.11
uv pip install --python .venv/bin/python -r requirements.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python -m pytest -q
.venv/bin/snakemake --snakefile workflow/Snakefile --cores 1 --dry-run
.venv/bin/snakemake --snakefile workflow/Snakefile --cores 1 --printshellcmds
.venv/bin/snakemake --snakefile workflow/Snakefile --cores 1 --dry-run
```

Acquisition requires HTTPS access to raw.githubusercontent.com on first run.
No private cache is required. Interrupted downloads are not accepted as inputs.
Corrupt cache files fail closed: explicitly remove the identified bad file
and rerun. Do not edit the pinned checksums to accommodate unverified bytes.

## Outputs (generated locally, ignored by Git)

- `data/cache/covid/`: six processed matrices and source metadata.
- `data/manifests/covid_all.csv`: all 114 source-order observations.
- `data/manifests/covid_primary.csv`: earliest available T1 then T2 per patient,
  85 rows before any later scientific exclusions; source order retained.
- `data/manifests/covid_sensitivity.csv`: all 114 observations with clinical ID
  as resampling group. These are **not train/test splits**; later resampling
  must keep each patient together.
- `results/manifests/{acquisition,audit,run_receipt,readiness}.json`: compact
  aggregate receipts and source/config/input/output hashes; no patient rows.

Metadata contain 85 patients and 29 paired patients. Only the local manifests
include necessary sample/patient IDs, timepoints and the unchanged source WOS
recoding, not demographics. Do not commit data, full metadata or model objects.
The receipt hashes the full local manifests without embedding their rows.

Expression matrices are processed cell-type means, **not raw counts or
count-aware pseudobulk**. No library-size normalization, count reconstruction,
biological interpretation or additional processing is implied. WOS 1.5 means
the original category **"1 or 2"**, not a fractional clinical severity score.
The source states selection optimized DIVAS correlation sensitivity; biological
and method-comparison generalizability may therefore be limited.

## Honest figure targets

`config/figures.yaml` registers exactly Fig1–Fig5, S1–S10, conditional S11/S12,
and TableS1–TableS8. Readiness records remain NOT_READY. For example:

```bash
.venv/bin/snakemake -s workflow/Snakefile --cores 1 results/figures/Fig1.pdf
```

This command intentionally fails with the missing scientific prerequisite path;
it never emits a placeholder PDF. Even if someone supplies a file at that path,
rendering still fails as NOT_IMPLEMENTED until a reviewed real renderer exists.
Conditional supplements are not approved for inclusion. This repository does
not meet the full paper reproducibility acceptance criterion yet.

## Licenses, attribution, future work

The MIT license covers only this independently authored foundation. The COVID
upstream states AGPL-3 for its code; a separate dataset license has **not** been
verified and is **not inferred** from that code license. No upstream code or
data are redistributed in this source tree/archive. Attribution and CLL audit:
`docs/PROVENANCE.md`. Future gates and adapter contract:
`docs/FUTURE_CONTRACT.md`. Reproduction, limitations and environment evidence:
`docs/REPRODUCIBILITY.md`.

## Pinned analysis package

The paper workflow runs against the [`rajiveplus`](https://github.com/mdmanurung/rajiveplus)
R package, pinned as a git submodule (`rajiveplus/`) at commit
`a94210ade9238553808c0d881824cc34c2c56f58` (stage-2, PR #9 head at pin time;
RJP-BASE v2, version 0.1.0). `config/package.yaml` records the pin and the
install procedure. Clone with `git clone --recurse-submodules` (or
`git submodule update --init`) to restore the exact package state; the
workflow builds and installs the package from that checkout before any
fitting. The pin tracks the master tip at tag v0.1.0-stage2.
