# SigBenchR — One-Page Lab Pitch

## Problem

Immunotherapy gene-signature papers often mix incompatible expression scales, silently recycle unstable gene symbols, and evaluate AUC on cohorts that still contain SD/NA labels used later for scoring—an invitation to leakage. Null permutation, multiplicity control, and small-n power are frequently missing, so “AUC = 0.72” is hard to interpret across anti-PD-1 RECIST and anti-CTLA-4 clinical-benefit definitions.

## Method

**SigBenchR** is a zero-Bioconductor, pure-R benchmarking layer. Datasets use `list(expr, sample_meta, info)` with optional GMT / `SummarizedExperiment` adapters. Five scorers (`mean`, `zscore`, `median`, `rank`, `ss_rank`) support signed signatures; 53 immuno-oncology aliases are rescued; coverage is labeled `ok` / `low_power` / `invalid`. SD/NA samples are dropped via `enc$kept` before scoring and permutation. Metrics include oriented DeLong ROC-AUC, stratified bootstrap PR-AUC, vectorized Mann–Whitney nulls with (k+1)/(B+1) and BH adjustment, REML + truncated HKSJ on logit-AUC, Hanley–McNeil power, biological-context axes (`correlate_reference_axes()`), and optional Cox PH / Kaplan–Meier survival (`evaluate_survival()`).

## Key Figure & Narrative

Read the normalization-labeled forest or heatmap **with** IFN-gamma axis correlation—not alone. Where IFNG6 tracks the IFN-gamma reference axis and retains AUC under anti-PD-1 RECIST cohorts, the biology is coherent; where correlation holds but discrimination collapses (or the reverse) under anti-CTLA-4 clinical-benefit labeling, the plot flags definition- and platform-driven breakdown rather than a single “winner” signature. Negative controls and coverage flags keep weak cells from looking definitive.

## Honest Limitations

Bundled demo matrices use a ~2,800-gene pool (full-matrix analyses use `read_user_dataset()` / `load_cohort_bundle()` with automatic `coverage_status`). Meta-analysis with `k <= 4` is exploratory. Cohort normalizations remain heterogeneous; SigBenchR never merges raw expression across platforms.

## Complementarity with TCGEx

TCGEx is a rich interactive exploration environment. SigBenchR is the lightweight, standalone, GMT- and SE-compatible benchmarking layer that stress-tests a candidate signature across cohorts with leakage controls, nulls, power, meta-summary, and optional survival—so exploratory discovery and confirmatory discrimination stay distinct.
