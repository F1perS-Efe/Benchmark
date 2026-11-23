# SigBenchR 0.4.5

* **Phase 1–7 Scientific, Statistical & Architectural Hardening:**
  * **Runtime & Reproducibility (Phase 1):** Unified deterministic `reproducibility_hash` (MD5) across standalone `export_manifest()`, `export_publication_bundle()`, and `export_report()`; locked Shiny reactive outputs and exports to the immutable `results()` snapshot; centralized `.resolve_min_genes()` so 2-gene `Cytolytic_Activity` (`GZMA`, `PRF1`) is scored consistently across all functions; added signed custom `genes_up`/`genes_down` R code generation in `generate_r_code()`.
  * **Cancer Immunology & Signature Fidelity (Phase 2):** Implemented continuous gene weight support (`weights`) across `score_signature()`, `run_benchmark()`, `analyze_gene_contributions()`, `null_permutation_test()`, and `generate_r_code()` using the published Ayers (2017) / Cristescu (2018) 18-gene `T_Cell_Inflamed_GEP` regression coefficients (including negatively weighted `CD276 = -0.023900` and `STAT1 = 0.250229`, with explicit disclaimers that RNA-seq evaluation does not claim analytical equivalence to the commercial NanoString nCounter assay); separated `Checkpoint_Axis` (6-gene classic PD-1/PD-L1/CTLA-4/LAG-3/TIGIT axis) from `B7_Immune_Exclusion` (`CD276`, `VTCN1`); removed pan-macrophage marker `CD68` and IFN-gamma target `IDO1` from `Myeloid_Suppression`; standardized 15 metadata fields across `get_builtin_signatures()`.
  * **Mouse-to-Human Ortholog & Background Integrity (Phase 3):** Expanded curated 1:1 MGI-HGNC ortholog dictionary to 160+ genes; added `ortholog_table`, `ortholog_policy = c("strict_1to1", "many_to_1")`, explicit 1:many ambiguity rejection, and provenance metrics (`attr(res, "ortholog_stats")`); preserved full transcriptome background in `attr(res, "transcriptome_bg_expr")` (`keep_unmapped_context = TRUE`) for `ss_rank` and `null_permutation_test()`.
  * **Statistical Engine, Null Permutation & Beyond Infiltration (Phase 4):** Eliminated CV test-set leakage in `test_beyond_infiltration()` via fold-inside gene scaling (`zscore`/`rank`) and predictor standardization; replaced pooled out-of-fold AUC with paired fold-level CV delta-AUC and uncertainty (`sd_delta_auc`, `delta_auc_ci_low`, `delta_auc_ci_high`); gated `adds_independent_value` on BH-adjusted `lrt_p_adj_bh < 0.05`; introduced sample-size weighted directional Stouffer Z synthesis (`combined_lrt_p`, `stouffer_z`, `direction_conflict`) alongside unsigned `fisher_omnibus_p`; added EPV (`epv`, `low_epv`) and VIF (`vif`, `high_collinearity`, `analysis_status`) guards; upgraded `null_permutation_test(match_expression = TRUE)` to 2D mean-quantile x SD-tertile matching with `background_type` classification; added `null_p_adj_bh_global` and `scope` to `adjust_null_p()`; separated `cox_p` (`p_adj_bh`), `logrank_p` (`logrank_p_adj_bh`), and unadjusted diagnostic `ph_p` in `run_survival_benchmark()`; hardened `pr_auc()` tie handling (`tie_fraction`, `small_positive_class`, `n_boot_valid`, `ci_unreliable`).
  * **Shiny UI/UX, Reporting & Resilience (Phases 5–8):** Integrated all Phase 4 diagnostic columns, human-readable explanations, collapsible statistical glossaries, paired CV delta-AUC forest plot, and rule-based warning banners into Shiny (`inst/shiny/app.R`), static HTML reports (`inst/rmd/benchmark_report.Rmd`), and `export_publication_bundle()`; added multi-cohort failure isolation, writable directory guards, and vectorized `zscore` matrix scoring; added methodology vignette `vignettes/SigBenchR-methodology.Rmd`.
* Documentation release: CI badge, overview and publication-ready
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
