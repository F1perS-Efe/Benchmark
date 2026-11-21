# SigBenchR 0.4.5

* Documentation release: CI badge, English overview and publication-ready
  Methods paragraph in `README.md`.
* Added CRAN-style `NEWS.md`, one-page lab pitch (`docs/LAB_PITCH.md`), and
  one-command HTML demo script (`inst/scripts/render_demo_report.R`).
* Ignored `docs/` and demo HTML artifacts in `.Rbuildignore` / `.gitignore`
  to keep CRAN check size and notes clean.
* P0 hygiene: ignore/untrack `.cursor/` and `sigbenchr_audit.zip`; align Shiny
  version badge with `DESCRIPTION` (0.4.5); fix stale README claim that gene
  alias rescue was absent.

# SigBenchR 0.4.4

* Added survival / clinical-benefit evaluation: `evaluate_survival()`,
  `run_survival_benchmark()`, and `combine_response_and_survival()`.
* Continuous Cox PH (1-SD scaled score) with HR, Wald p, and C-index; median
  split Kaplan–Meier log-rank with group median survival.
* Optional `os_time`/`os_event` and `pfs_time`/`pfs_event` columns validated in
  `validate_dataset()` (backward compatible when absent).
* `export_report()` accepts optional `surv_res`; report renders a Clinical
  Benefit section when survival results are present.
* Added Recommended package `survival` to `Suggests` (loaded via
  `requireNamespace()`).

# SigBenchR 0.4.3

* Biological context axes via `correlate_reference_axes()` (IFN-gamma,
  cytotoxicity, T-cell infiltration).
* Hanley–McNeil power analysis (`estimate_auc_power()`) with per-row
  `power_auc70` in `run_benchmark()`.
* Normalization labels and stratified reporting helpers for forest/heatmap
  interpretation across heterogeneous ICI cohorts.

# SigBenchR 0.4.2

* Gene alias rescue for 53 immuno-oncology symbols; `coverage_status`
  (`ok` / `low_power` / `invalid`) with clearer low-coverage messaging.
* Expanded built-in literature ICI panel (8+) plus negative controls via
  `get_builtin_signatures()`.
* GMT adapter (`read_gmt()`) and optional `SummarizedExperiment` coercion
  without hard Bioconductor Imports.

# SigBenchR 0.4.1

* Fixed ROC SD leakage (score/encode only on `enc$kept` samples).
* Shiny init/report synchronization fixes; signed-null edge cases; README
  honesty updates for bundled gene-pool limits.

# SigBenchR 0.4.0

* Logit-AUC REML + truncated HKSJ random-effects meta-analysis in
  `summarize_auc()` / forest plots.
* SD-handling spectrum sensitivity via `check_sd_sensitivity()`.
* Benchmark vignette for reproducible ICI signature evaluation.

# SigBenchR 0.3.1

* Redesigned Shiny layout: collapsible controls, KPI strip, contextual tabs.

# SigBenchR 0.3.0

* Cohort-independent `ss_rank` scoring (pure base R).
* Signed `+`/`-` signatures (`genes_up` / `genes_down` or named numeric genes).
* Subgroup stratification (`stratify_by`) for treatment / endpoint slices.

# SigBenchR 0.2.0

* Local user dataset contract (`read_user_dataset()`).
* Lossless cohort bundle loader (`load_cohort_bundle()`).
* Detailed gene coverage reporting for arbitrary signatures.

# SigBenchR 0.1.1

* Data provenance documentation and cohort integrity tests.
* GitHub Actions `R-CMD-check` workflow.

# SigBenchR 0.1.0

* Initial CRAN-oriented MVP: `list(expr, sample_meta, info)` contract,
  multi-method scoring, ROC-AUC / PR-AUC, null permutation, weighted
  summaries, Shiny UI, and curated melanoma ICI demo cohorts.
