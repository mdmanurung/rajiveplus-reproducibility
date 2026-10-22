# Reproduction and deviations

Use the README commands from the repository root. Python 3.11.13 was observed.
`requirements.txt` pins only actually observed direct workflow dependencies.
`environment.yml` is an installation recipe, not a solved lockfile. No
`renv.lock`, fabricated conda lock, R-method environment or historical
reproduction claim is provided. The accompanying external observed package
inventory records the complete installed Python environment used for testing.

Snakemake uses the **local executor**, `--cores 1`; SLURM is unavailable.
No scheduler job IDs, inference seeds, fits or benchmark result rows exist.
Use single-thread environment variables shown in README. Run receipts record
actual Python/dependency/platform versions, thread settings, Git identity when
available, source/config/input/output hashes, audit time and process RSS.
RSS is the audit process high-water mark on Linux, not total DAG peak memory
or a method-performance measurement. Acquisition time may include cache
verification rather than downloading; command logs distinguish clean runs.

Forcing a complete smoke check after a config/source change:

```bash
.venv/bin/snakemake -s workflow/Snakefile --cores 1 --forceall
```

The source archive lacks `.git`; there the Git fields truthfully report
`UNAVAILABLE_SOURCE_ARCHIVE`. The external handoff identifies the exact
committed source and archive checksum. Config/input hashes remain available.
To preserve local history, use the accompanying git bundle instead.

The receipt's source hashes cover executable workflow source and requirements;
all five configuration files have separate hashes. Snakemake's dependencies
invalidate downstream receipts when these files change. Receipts include
timestamps/timings so they are not promised bitwise identical across machines;
input/config/manifest identities are deterministic.

Unit tests use tiny invented **plumbing-only** fixtures and never create
scientific results. Scientific readiness is written separately from plumbing
PASS: all methods, comparisons, simulations, CLL, calibration and clinical
analysis remain NOT_RUN. Five main scientific figures and all supplements/
tables cannot yet be reproduced; requests fail rather than generating fake
graphics. No manuscript or empty comparator wrapper is shipped.

Data, local manifests, results and logs are ignored and excluded from source
archives. Retain a separate local cache if needed, but do not redistribute
clinical rows without checking terms. The compact final audit/run receipts
contain hashes and aggregate counts, not per-patient rows.
