# Future work is NOT_IMPLEMENTED / NOT_RUN

This is a pilot plumbing milestone, not a confirmatory protocol or complete
scientific repository. The full publication specification remains the target.
No fake historical manifests or substitute benchmark values are introduced.

Future method adapters must accept named sample-by-feature matrices, ordered
sample/patient IDs, observed masks, scenario config and an isolated seed. They
must return a terminal status (success, failed, nonconverged, timeout, OOM,
unsupported), estimator-specific fit and only supported subspaces,
reconstructions and weights. Records must link method source/version,
input/config hashes, pilot/confirmatory class, seed, elapsed time, peak RSS,
exit code and complete logs. Unsupported outputs must not be mapped to fake
common metrics. No adapter is implemented in this milestone.

Seed policy before future simulations: freeze scenario IDs, seed derivation,
replicate allocation, certified truths and held-out confirmatory manifests.
No numeric seed is listed here because no random operation was executed.
Projector loss, exact rank, over-selection and Type-I error require precise
metric definitions and uncertainty before comparisons.

Do not weaken the existing future gates:

- geometry Gram error <=1e-10; Individual leakage <=1e-8;
- rank null: 200 datasets, Type-I <=0.10;
- rank one: 100 datasets, exact recovery >=0.80;
- rank two: 100 datasets, exact recovery >=0.80, over-selection <=0.10;
- original candidate projector-loss degradation margin +0.02, contingent
  on recovering its original metric definition (not redefined here).

Wedin/geometry/rank instruments, source/environment pins, calibrated inference
and generator truth must be resolved before promotion. Missing historical
gate manifests remain missing. New primary COVID analysis uses earliest
available sample per patient; sensitivity keeps all 114 with grouped
resampling. Never split a patient across train/test. Freeze endpoints,
covariates, multiplicity and ambiguous WOS handling before inference.
Selection optimized for DIVAS sensitivity may restrict biological findings
to exploratory interpretation. Neither application assumes a winning method.

Next workflow expansion must add real implementations and explicit scientific
prerequisites, not empty scripts or dummy result files. The 25 registered
figure/table items are identifiers only, not claims of completed designs.
