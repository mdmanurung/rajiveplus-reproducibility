# Input provenance, semantics and independent authorship

## COVID processed-input track

Source: ByronSyun/DIVAS_COVID19_CaseStudy, commit
`12cae2bd385919cb2ffb48abb07813be9c440fc3`, directory
`scRNA_celltyist_analysis/DIVAS_run/divas_input_T1T2combined`.
The typo `celltyist` is intentional and matches upstream.
`config/datasets.yaml` contains exact SHA256s measured from the pinned bytes,
then independently verified by the executable workflow.

Pinned documentation:

- `https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy/blob/12cae2bd385919cb2ffb48abb07813be9c440fc3/scRNA_celltyist_analysis/README.md`
  describes processed matrix dimensions, filtering, timepoint balance, and
  selection optimized for DIVAS correlation sensitivity. It does not imply
  independence of all 114 observations.
- `https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy/blob/12cae2bd385919cb2ffb48abb07813be9c440fc3/scRNA_celltyist_analysis/T1T2_preparation/combine_t1_t2_data.R`
  aggregates expression with `mean(axis=0)` from processed AnnData; means are
  not raw counts. No raw-cell reconstruction was executed here.
- Source metadata `final_sample_info.csv` supplies source-order sample IDs,
  clinical IDs, T1/T2 and WOS recodings. The original WOS `"1 or 2"` recoding
  to 1.5 is retained as an ambiguous category, per the approved plan; it must
  not be interpreted as a measured fractional clinical score.
- Upstream `README.md` and `LICENSE` state AGPL-3 for the source project.
  Separate data/clinical-metadata terms remain unverified. Public availability
  of IDs is not a reason to redistribute unnecessary patient-level metadata.

Six matrices must align *exactly* to metadata sample order; set equality alone
is insufficient. No silent reordering, transposition, imputation or zero-filling.
Checks cover bytes, row widths/counts, duplicate/empty IDs/features, finite
numeric values, patient/timepoint uniqueness, patient count, pair count and
T1/T2 counts. Noninteger expression is an additional semantic consistency
check, not by itself proof of a particular upstream normalization.

Manifests are selection indexes only. All 114 source rows are preserved in
the all/sensitivity indexes; no demographics are exported to these indexes.
Primary selects the earliest **available** timepoint, not necessarily baseline
for everyone and not the first row encountered. Clinical ID is the grouping
unit in both arms. No folds or inferential endpoints have been chosen.

## CLL vignette audit: NOT_RUN acquisition

Reviewed package source:
`mdmanurung/rajiveplus@8ceabeedea6501a15b10158e40ad8cebd904ec11`,
`vignettes/cll_application.Rmd`. No CLL package/data download or fit executed.
The vignette cites Dietrich et al. (2018), PMID 29227291 and
BloodCancerMultiOmics2017, with the original MOFA/MOFA2 tutorial as context.
Its input objects are `dds`, `methData`, `lpdAll`, and `patmeta`. Its code:

1. Maps RNA columns using `colData(dds)[["PatID"]]` to patient H-IDs.
2. Transposes `log2(counts + 1)` and selects 500 high-variance RNA features.
3. Transposes methylation, selecting 500 high-variance CpGs.
4. Filters drug features to `type == "viab"`, transposes and chooses up to
   448 variable features.
5. Takes the three-view intersection, then restricts `patmeta$Diagnosis`
   to CLL, excluding missing diagnosis; centers/scales each block.
6. Its scaling helper replaces NaN with zero; this must be audited explicitly
   for zero-variance features and missingness rather than copied silently.

The approximate ~121 overlap is a vignette claim, **not verified here**.
Before acquisition/analysis: pin the Bioconductor package/data version and
accession, verify data-specific terms independently of software licensing,
validate uniqueness and completeness of all ID mappings, freeze inclusion/
exclusion and preprocessing, preserve missingness masks, and audit whether
feature selection must be training-only for a prediction endpoint. Binary
mutations are not a Gaussian input view in that vignette. This stage does not
adopt its clinical models or reproduce its analysis.

## Software attribution

All stage-0 workflow code, tests and prose in this repository were newly
authored for the foundation, not copied from the package or upstream
COVID/CLL source. Upstream algorithms remain external, pinned in
`config/methods.yaml`; their independent licenses and citations remain
applicable. There has been no benchmark-script migration or relicensing.
