# Stage 1 paired geometry and Individual-handling ablation protocol

**Protocol status:** `DESIGN_AND_SMOKE_VALIDATED_NO_BENCHMARK`
**Scientific readiness:** `NOT_READY`
**Calibration, simulation study, method promotion, figures, and manuscript
claims:** `NOT_RUN`

This is a new prospective Stage 1 design, not a reconstruction of a historical
preregistration. It defines a nonconfounded paired \(2\times2\times2\)
ablation and a fixed-data smoke-test contract. It does not authorize method
fits, automatic-rank calibration, scientific comparisons, or figure
generation. The machine-readable source of truth is
[`config/stage1_ablation.json`](../config/stage1_ablation.json).

## 1. Estimand and eight paired arms

Every dataset is evaluated in all eight arms. Dataset matrix bytes, block
order, preprocessing, initial block signal ranks, requested Joint rank,
contamination mask, cached robust block fits, and applicable randomization
identifiers are held fixed. Factor codes are:

| factor | level 0 | level 1 |
|---|---|---|
| **R — representation** | raw block scores and raw final Joint scores | rank-revealing orthonormal block scores and final Joint scores |
| **S — stack aggregation** | robust stacked-score factorization | classical SVD/eigenanalysis of the same represented scores |
| **I — Individual handling** | original selected robust Individual reconstruction | that same reconstruction projected out of the fitted Joint span and algebraically refactored |

The arms are `R0S0I0`, `R0S0I1`, `R0S1I0`, `R0S1I1`, `R1S0I0`,
`R1S0I1`, `R1S1I0`, and `R1S1I1`. `R0S0I0` is the frozen baseline and
must reproduce baseline correctness fixtures. A classical **block** fit is
not a primary factor and must not be substituted into any S=1 arm.

### R: block and final representations are both explicit

R=0 passes each cached robust block-score matrix to aggregation unchanged.
After aggregation it preserves raw Joint-score coordinates unchanged,
including their order and signs. It does not silently
QR-normalize either representation.

R=1 applies a deterministic rank-revealing thin orthonormalization to each
block-score matrix before aggregation and to the selected final Joint scores
after aggregation. Numerical rank, tolerance, pivots/order, and sign
orientation must be recorded. A deficient requested representation is
reported as rank-deficient; columns are not silently invented.

For both R levels, this representation-specific **final-score
canonicalization precedes the common l2 identifiability filter**. The filter
implementation, tolerance, and ordering are identical in every arm. Outputs
must include canonical pre-filter scores/rank, retained and dropped direction
IDs, post-filter scores/rank, and pre/post Gram diagnostics. Thus R changes
the declared block and final representations, while the l2 rule remains a
common downstream operation.

### S: aggregation does not change block fitting

`.ablation_prepare` performs the robust block fits once and caches them. S=0
uses the current robust stacked-score factorization. S=1 uses ordinary
stacked SVD (or the algebraically equivalent eigenanalysis) on the **same
represented cached scores**. It never computes classical block SVDs.

For R=0 the classical Gram is an amplitude-weighted raw-score Gram and is not
a projector sum. Only R=1 admits the projector-sum interpretation. This
distinction is retained in output metadata rather than being hidden by an
unconditional orthonormalization.

### I: one selected robust reconstruction, then an algebraic constraint

The robust Individual reconstruction and its selected rank are obtained once
after each fixed R/S geometry and its identifiability filter.
I=0 returns that reconstruction.
I=1 computes

\[
  \widetilde I_k=(I-Q_JQ_J^\mathsf{T})\widehat I_k
\]

using an orthonormal basis Q_J of the fitted Joint **span**, even for raw
scores, and then algebraically refactors
\(\widetilde I_k\). It does not robust-refit. The pre-constraint selected rank
is frozen within the I pair. Each block records that rank, post-projection
numerical rank, and rank drop. Geometry is immutable during finishing, so an I
pair shares an identical fitted Joint geometry and an identical preselected
Individual reconstruction.

## 2. Planned internal interface and ordered execution

The package-side experimental helpers are implemented as internal-only:

```text
.ablation_prepare(blocks, initial_ranks, seed)
.ablation_geometry(cache, joint_rank,
                   representation='raw'/'orthonormal',
                   aggregation='robust'/'classical', seed)
.ablation_select_individual(cache, geometry, seed)
.ablation_finish(cache, geometry,
                 individual='original'/'constrained', seed, selected)
```

They must never route into the public estimator default. The order is:

1. validate fixed inputs and certified truth;
2. fit/cache robust blocks and thresholds once;
3. apply the R block representation;
4. apply S aggregation;
5. canonicalize the R-specific final Joint scores;
6. emit pre-l2 scores, ranks, and diagnostics;
7. apply the common l2 identifiability filter;
8. emit post-l2 scores, ranks, retained/dropped directions, and diagnostics;
9. select Individual once per R/S geometry, then finish both I levels from
   the shared selection and frozen geometry; and
10. emit terminal status and all rank-drop/failure records.

Seeds are isolated by stage. The executed smoke uses explicit data/block/
aggregation/Individual seeds 1729/8101/8102/8103. A SHA-256 label derivation is
specified for future runs. Individual-selection labels depend on R/S but
never I. Null-permutation labels are independent of arm so
the underlying permutations can be paired, while block labels are distinct
so permutations are independently generated within each block.

## 3. Rank layers and blocked alternatives

### Authorized first layer: fixed rank

The first layer uses certified synthetic truth for the initial block ranks and
target Joint rank. It reports geometry before and after the l2 filter. This is
the only rank layer authorized for the planned smoke target.

### Automatic rank: BLOCKED

Automatic rank selection is `BLOCKED_PENDING_MATCHED_NULL`. Its candidate null
independently permutes rows within each block, computes each arm's own matched
statistic, and reuses the same block-specific underlying permutations across
all arms. Permutation generation has an isolated stage seed and is not
interleaved with fitting RNG. Draw count, tail rule, and calibration remain
unfrozen; therefore this layer cannot run.

An arm that has no effective finite threshold returns
`BLOCKED_NO_EFFECTIVE_FINITE_THRESHOLD`, never a successful fit based on a
negative-infinity or missing cutoff.

The Wedin candidate is a separate named alternative, not a fallback, maximum,
or minimum combined with the permutation candidate. Cross-reference the
companion package work described by approved Stage 1 plan section 2,
“Source-grounded Wedin derivation.” This paper repository makes no theorem or
bibliographic claim and includes no inherited source bundle. The companion
`stage1-wedin-derivation.md` resolves the classical algebra and pinned-source
audit, not robust transfer or coverage. The alternative remains
`BLOCKED_UNRESOLVED_ROBUST_THEOREM_AND_COVERAGE`.

## 4. Future synthetic scenarios (specified, not run)

All core scenarios use three blocks with \(n=60\) and
\(p=(30,25,20)\). Use the existing centered generator
`rajiveplus_identifiable_v1` (`.simulate_identifiable_data`) with `snr=Inf`,
Joint strengths `sqrt(n)*joint_amplitudes` and Individual strengths
`sqrt(n)*1.25`. It transforms a Gaussian score pool by
\(H=I-\mathbf1\mathbf1^\mathsf{T}/n\), then thin-QR orthonormalizes it.
Joint scores are shared; disjoint pool columns give mutually orthogonal
block-specific Individual scores in the Joint complement. Combined
Joint/Individual feature loadings are orthonormal within each block.
Component singular values are \(\sqrt n\) times the amplitudes
below, and centered noise is \(HZ_k\), with iid pre-transform
\(Z_{k,ij}\sim N(0,1)\). This separately seeded unit-noise scenario adapter
is specified, not implemented or run in Stage 1; it is not the SNR-scaled smoke.

| scenario | Joint rank | Joint amplitudes | Individual ranks | Individual amplitude |
|---|---:|---|---|---:|
| `joint_rank0` | 0 | `[]` | `(1,1,1)` | 1.25 |
| `joint_rank1` | 1 | `(2.0)` | `(1,1,1)` | 1.25 |
| `joint_rank2` | 2 | `(2.0,1.5)` | `(1,1,1)` | 1.25 |
| `joint_rank2_weak_last` | 2 | `(2.0,0.35)` | `(1,1,1)` | 1.25 |

These amplitudes define fixtures, not calibrated effect sizes. Initial block
ranks are Joint rank plus one. No failed draw is silently resampled. After the
centering transform, an independent matrix-based truth certifier must verify
component/noiseless ranks, Joint/Individual orthogonality, zero all-block
Individual intersection, the exact all-block Joint intersection, finite
values, and reconstruction. A failed certificate is a terminal invalid-truth
record.

Each latent draw has paired contamination conditions:

* `clean`: no altered entries.
* `cell_5pct_plus10_robust_scale`: in each block select exactly 5% of
  cells—90, 75, and 60 cells—uniformly without replacement. Add
  \(+10s_k\), where
  \(s_k=\operatorname{median}|E_k-\operatorname{median}(E_k)|/
  0.6744897501960817\), computed before contamination. The mask and signs are
  fixed across all eight arms.
* `row_10pct_plus10_robust_scale` is an optional development-only sensitivity:
  choose exactly six sample rows, shared across blocks, and add \(+10s_k\) to
  every cell in those rows. It has no confirmation allocation or gate.

Contamination is added after centered latent truth has been constructed and
certified; adapters do not recenter it. This keeps both the post-transform
truth and clean-entry mask unambiguous.

## 5. Finite seed manifest and pairing

All listed pools are finite explicit integer arrays in the JSON manifest.
Only the separate smoke seed was consumed; no study pool was used.

* smoke fixture: seed `1729`, separate from study pools, n=36, p=(12,10,8),
  Joint strengths (4,3), Individual strength 2, SNR=4;
* development: `1001` through `1020` (20 seeds), future and not run;
* untouched confirmation null: `20001` through `20200` (200 seeds);
* untouched confirmation rank one: `21001` through `21100` (100 seeds);
* untouched confirmation rank two: `22001` through `22100` (100 seeds).

A seed identifies one latent dataset and all of its paired arms/declared
contamination versions. Confirmation pools are locked against development,
tuning, threshold choice, or smoke use. The weak-last-signal and optional-row
conditions are development-only unless a later protocol revision freezes new
confirmation pools before inspection.

## 6. Endpoints, contrasts, uncertainty, and failures

The historical primary rank-two endpoint is retained exactly as

\[
 L_{\mathrm{hist},r=2} =
 \lVert P_{\mathrm{est}}-P_{\mathrm{true}}\rVert_F/\sqrt{4}.
\]

The candidate degradation margin `+0.02` applies **only** to that exact
rank-two endpoint and normalization. It is not transferred to rank zero,
rank one, weak-signal, unnormalized, reconstruction, or robust-loss
endpoints.

P_est is always the **orthogonal projector onto the numerical fitted span**,
recomputed from a rank-revealing basis for every arm; it is never raw U U^T
for nonorthogonal robust scores. Because I is applied after immutable Joint
geometry, I contrasts on Joint projector/rank endpoints are structural zeros,
not evidence that the Individual constraint is ineffective. I is evaluated
using Individual leakage and reconstruction/robust-loss endpoints. Projection
is an algebraic constraint, not a constrained robust optimum, and can worsen
other losses.

Other recorded endpoints are actual pre/post ranks; exact, under-, and
over-selection; block/final Gram errors; algebraic Individual leakage into the
fitted Joint span; leakage into the true Joint span; normalized reconstruction
loss; a fully specified Huber robust loss; and clean-entry reconstruction
error when a contamination mask is known. The Huber scale is the
precontamination noise MAD divided by 0.6744897501960817, shared across arms;
zero scale is unsupported. The algebraic constrained-arm
leakage gate and true-span scientific leakage are distinct.

For endpoint value \(y_d(r,s,i)\), code levels as \(x=-1,+1\). For every
nonempty subset \(A\) of `{R,S,I}`, the dataset-level factorial effect is

\[
  e_{d,A}=2^{-2}\sum_{r,s,i}
  \left(\prod_{a\in A}x_a\right)y_d(r,s,i).
\]

This yields the three paired main effects and `RS`, `RI`, `SI`, and `RSI`
interactions. Uncertainty resamples datasets while retaining each eight-arm
cluster. The planned interval is a 95% percentile cluster bootstrap with
10,000 resamples; Holm adjustment is within each declared endpoint family
across the seven contrasts. These are future analysis settings, not executed
in this milestone.

For the exact rank-two primary loss, separately report each candidate-minus-
baseline paired mean difference. Nondegradation requires its one-sided
Bonferroni upper percentile-bootstrap bound (family alpha 0.05 across seven
candidate arms and two core contamination conditions) to be below +0.02.
No superiority is inferred from passing this margin. Descriptive factorial
intervals are nominal 95%; any factorial inferential p-values require a later
frozen test construction before confirmation data are accessed. Holm is the
declared adjustment, not a claim that percentile intervals themselves are
Holm-adjusted. All inference remains blocked until automatic-rank calibration
and the separately approved experiment implementation are complete.

Every attempted dataset-arm cell contributes to a terminal-status table:
`success`, `failed`, `nonconverged`, `rank_deficient`, `invalid_truth`,
`unsupported`, `timeout`, or `oom`. Failure indicators are analyzed over all
arms. Numeric contrasts state both the all-attempted denominator and complete
eight-arm count; incomplete octets and missing-arm patterns are reported, not
silently discarded. With rank-two truth and estimated-rank capacity three,
the historical projector loss is bounded by sqrt(5)/2, not one (one requires
equal estimated/true rank two). The failure sensitivity maps failed rank-two
arms to sqrt(5)/2. Unbounded endpoints are
not imputed and cannot support superiority claims when informative failures
remain.

## 7. Gates and interpretation

Existing future gates are unchanged:

* block and final orthonormal Gram error at most `1e-10`;
* constrained Individual leakage into the fitted Joint span at most `1e-8`;
* null: 200 datasets and Type-I error at most `0.10`;
* rank one: 100 datasets and exact recovery at least `0.80`;
* rank two: 100 datasets, exact recovery at least `0.80`, and over-selection
  at most `0.10`.

The rank gates cannot be evaluated while automatic rank is blocked. Gram and
leakage gates are mechanical correctness checks, not evidence that an arm is
scientifically superior.

## 8. Executable boundary

The fixed-data eight-arm smoke orchestrator is implemented at
`scripts/run_stage1_smoke.R` and passed on one certified centered synthetic
dataset. Raw/robust/original matches public baseline scores and numerical
component matrices; public dimension names are checked separately. The
scientific study, unit-noise adapter and automatic-rank layer have not run.
Smoke execution validates arm
reachability, cache sharing, ordering, baseline equivalence, rank logging, and
terminal schemas only. It may not estimate performance, calibrate a threshold,
choose a winning arm, alter public defaults, or populate scientific figures.

All manuscript software and figures remain `NOT_READY`. No inherited source
files or source bundles are included.
