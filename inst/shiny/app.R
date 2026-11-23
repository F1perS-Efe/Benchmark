# SigBenchR Shiny - tek dosya

# Global: oturum basina yeniden yuklemeyi onlemek icin server disinda
data(datasets_ici, package = "SigBenchR", envir = environment())
data(signatures_demo, package = "SigBenchR", envir = environment())

.builtin_sigs_ui <- SigBenchR::get_builtin_signatures(include_controls = TRUE)
.builtin_sigs_pool <- c(
  .builtin_sigs_ui,
  signatures_demo[setdiff(names(signatures_demo), names(.builtin_sigs_ui))]
)
.sigbenchr_version <- tryCatch(
  as.character(utils::packageVersion("SigBenchR")),
  error = function(e) "0.4.5"
)

# Fail-closed: yalniz run_app(allow_upload = TRUE) ile acilir
.sigbenchr_allow_upload <- isTRUE(
  shiny::getShinyOption("sigbenchr_allow_upload", default = FALSE)
)

.validate_custom_signature <- function(genes_text = NULL,
                                      genes_up = NULL,
                                      genes_down = NULL,
                                      custom_name = "Custom",
                                      existing_names = character(0)) {
  nm <- trimws(as.character(custom_name)[1])
  if (is.na(nm) || !nzchar(nm)) {
    nm <- "Custom"
  }
  if (length(existing_names) > 0L && tolower(nm) %in% tolower(as.character(existing_names))) {
    stop(
      "Custom signature name '", nm,
      "' already exists in selected builtin signatures (name collision).",
      call. = FALSE
    )
  }

  has_updown <- !is.null(genes_up) || !is.null(genes_down)
  if (has_updown) {
    gu <- if (is.null(genes_up)) {
      character(0)
    } else {
      SigBenchR::parse_genes(paste(as.character(genes_up), collapse = " "))
    }
    gd <- if (is.null(genes_down)) {
      character(0)
    } else {
      SigBenchR::parse_genes(paste(as.character(genes_down), collapse = " "))
    }
    if (length(gu) == 0L && length(gd) == 0L) {
      stop("Enter at least one valid gene for the custom signature.", call. = FALSE)
    }
    ov <- intersect(toupper(gu), toupper(gd))
    if (length(ov) > 0L) {
      stop(
        "genes_up and genes_down overlap: ",
        paste(ov, collapse = ", "),
        call. = FALSE
      )
    }
    return(list(name = nm, genes = NULL, genes_up = gu, genes_down = gd))
  }

  genes <- SigBenchR::parse_genes(if (is.null(genes_text)) "" else genes_text)
  if (length(genes) == 0L) {
    stop("Enter at least one valid gene for the custom signature.", call. = FALSE)
  }
  list(name = nm, genes = genes, genes_up = NULL, genes_down = NULL)
}

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$script(shiny::HTML(
      "function sigbenchrCopyRcode() {
         var el = document.getElementById('r_code_text');
         if (!el) return;
         var txt = el.value || '';
         if (navigator.clipboard && navigator.clipboard.writeText) {
           navigator.clipboard.writeText(txt);
         } else {
           el.focus();
           el.select();
           document.execCommand('copy');
         }
       }"
    )),
    shiny::tags$style(shiny::HTML("
      :root {
        --sb-bg: #FFFFFF;
        --sb-panel: #F8FAFC;
        --sb-surface: #FFFFFF;
        --sb-border: #E2E8F0;
        --sb-accent: #0F3D5E;
        --sb-accent-hover: #0A2F49;
        --sb-text: #1E293B;
        --sb-muted: #475569;
        --sb-amber-bg: #FFFBEB;
        --sb-amber-border: #F59E0B;
        --sb-amber-text: #92400E;
        --sb-radius: 6px;
      }
      body, .container-fluid {
        background: var(--sb-bg) !important;
        color: var(--sb-text);
      }
      .sb-header {
        display: flex;
        align-items: flex-start;
        justify-content: space-between;
        gap: 16px;
        padding: 14px 18px;
        margin: 0 0 12px 0;
        background: var(--sb-surface);
        border: 1px solid var(--sb-border);
        border-radius: var(--sb-radius);
      }
      .sb-brand-title {
        margin: 0;
        font-size: 1.55rem;
        font-weight: 700;
        color: var(--sb-accent);
        letter-spacing: -0.02em;
      }
      .sb-brand-sub {
        margin: 4px 0 0 0;
        font-size: 0.92rem;
        color: var(--sb-muted);
      }
      .sb-header-meta {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        align-items: center;
        justify-content: flex-end;
      }
      .sb-badge {
        display: inline-block;
        padding: 4px 10px;
        border-radius: 4px;
        font-size: 0.78rem;
        font-weight: 600;
        border: 1px solid var(--sb-border);
        background: var(--sb-panel);
        color: var(--sb-muted);
        white-space: nowrap;
      }
      .sb-badge-warn {
        background: var(--sb-panel);
        border-color: var(--sb-border);
        color: var(--sb-accent);
      }
      .sb-stale {
        margin: 0 0 12px 0;
        padding: 10px 14px;
        background: var(--sb-amber-bg);
        border: 1px solid var(--sb-amber-border);
        border-radius: var(--sb-radius);
        color: var(--sb-amber-text);
        font-size: 0.92rem;
      }
      .sb-kpi-row {
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: 10px;
        margin: 0 0 14px 0;
      }
      @media (max-width: 992px) {
        .sb-kpi-row { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      }
      @media (max-width: 600px) {
        .sb-kpi-row { grid-template-columns: 1fr; }
        .sb-header { flex-direction: column; }
      }
      .sb-kpi {
        background: var(--sb-surface);
        border: 1px solid var(--sb-border);
        border-radius: var(--sb-radius);
        padding: 12px 14px;
      }
      .sb-kpi-label {
        font-size: 0.75rem;
        text-transform: uppercase;
        letter-spacing: 0.04em;
        color: var(--sb-muted);
        margin: 0 0 6px 0;
      }
      .sb-kpi-value {
        font-size: 1.05rem;
        font-weight: 650;
        color: var(--sb-text);
        margin: 0;
        word-break: break-word;
      }
      .sb-welcome {
        background: var(--sb-surface);
        border: 1px solid var(--sb-border);
        border-radius: var(--sb-radius);
        padding: 16px 18px;
        margin: 0 0 14px 0;
      }
      .sb-welcome h4 { margin-top: 0; color: var(--sb-accent); }
      .sb-welcome ol { margin-bottom: 0; color: var(--sb-muted); }
      .well, .sidebarPanel, .mainPanel > .tabbable {
        background: var(--sb-surface) !important;
        border: 1px solid var(--sb-border) !important;
        border-radius: var(--sb-radius) !important;
        box-shadow: none !important;
      }
      .sidebarPanel {
        padding: 14px !important;
      }
      .sb-actions {
        display: flex;
        flex-direction: column;
        gap: 8px;
        margin-bottom: 12px;
        padding-bottom: 12px;
        border-bottom: 1px solid var(--sb-border);
      }
      .sb-actions .btn-primary,
      .sb-actions .btn-primary:focus {
        background-color: var(--sb-accent) !important;
        border-color: var(--sb-accent) !important;
      }
      .sb-actions .btn-primary:hover {
        background-color: var(--sb-accent-hover) !important;
        border-color: var(--sb-accent-hover) !important;
      }
      .sb-actions .shiny-download-link {
        text-align: center;
      }
      details.sb-acc {
        background: var(--sb-panel);
        border: 1px solid var(--sb-border);
        border-radius: var(--sb-radius);
        padding: 8px 10px;
        margin-bottom: 10px;
      }
      details.sb-acc > summary {
        cursor: pointer;
        font-weight: 650;
        color: var(--sb-text);
        outline: none;
        padding: 4px 2px;
      }
      details.sb-acc[open] > summary {
        margin-bottom: 8px;
        border-bottom: 1px solid var(--sb-border);
        padding-bottom: 8px;
      }
      .sb-hint {
        font-size: 0.82rem;
        color: var(--sb-muted);
        margin: 0 0 8px 0;
      }
      .sb-note {
        font-size: 0.84rem;
        color: var(--sb-muted);
        margin: 6px 0 0 0;
      }
      .sb-panel-block {
        background: var(--sb-surface);
        border: 1px solid var(--sb-border);
        border-radius: var(--sb-radius);
        padding: 12px 14px;
        margin-bottom: 12px;
      }
      .sb-inline-controls {
        display: flex;
        flex-wrap: wrap;
        gap: 12px;
        align-items: flex-end;
        margin-bottom: 10px;
      }
      .sb-inline-controls .form-group {
        margin-bottom: 0 !important;
        min-width: 160px;
      }
      .nav-tabs > li > a { color: var(--sb-muted); }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover,
      .nav-tabs > li.active > a:focus {
        color: var(--sb-accent) !important;
        border-top-color: var(--sb-accent) !important;
        font-weight: 600;
      }
      h4.sb-section {
        margin-top: 4px;
        margin-bottom: 10px;
        font-size: 1.05rem;
        color: var(--sb-text);
      }
      input[type='checkbox'], input[type='radio'] {
        accent-color: #0F3D5E;
      }
      .sb-alert-strip {
        display: flex;
        flex-direction: column;
        gap: 6px;
        margin: 0 0 12px 0;
      }
      .sb-alert-item {
        padding: 8px 12px;
        border-radius: var(--sb-radius);
        border: 1px solid #FDE68A;
        border-left: 4px solid #F59E0B;
        background: #FFFBEB;
        color: #78350F;
        font-size: 0.88rem;
      }
      .sb-alert-item.sb-alert-info {
        border-color: #BFDBFE;
        border-left-color: #3B82F6;
        background: #EFF6FF;
        color: #1E3A8A;
      }
      .sb-alert-item.sb-alert-danger {
        border-color: #FECACA;
        border-left-color: #EF4444;
        background: #FEF2F2;
        color: #991B1B;
      }
      .sb-alert-item.sb-alert-success {
        border-color: #A7F3D0;
        border-left-color: #10B981;
        background: #ECFDF5;
        color: #065F46;
      }
      .sb-stage-bar {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        align-items: center;
        justify-content: space-between;
        padding: 8px 12px;
        margin: 0 0 10px 0;
        background: var(--sb-panel);
        border: 1px solid var(--sb-border);
        border-radius: var(--sb-radius);
        font-size: 0.84rem;
        color: var(--sb-muted);
      }
      .sb-tier-pill {
        display: inline-block;
        padding: 2px 8px;
        border-radius: 4px;
        font-size: 0.8rem;
        font-weight: 600;
        border: 1px solid var(--sb-border);
        background: var(--sb-panel);
      }
      .sb-tier-high { background: #ECFDF5; border-color: #10B981; color: #065F46; }
      .sb-tier-med  { background: #EFF6FF; border-color: #3B82F6; color: #1E40AF; }
      .sb-tier-low  { background: #F8FAFC; border-color: #94A3B8; color: #334155; }
    "))
  ),
  shiny::div(
    class = "sb-header",
    shiny::div(
      shiny::tags$h1(class = "sb-brand-title", "SigBenchR"),
      shiny::tags$p(
        class = "sb-brand-sub",
        "Immunotherapy Gene Signature Benchmarking"
      )
    ),
    shiny::div(
      class = "sb-header-meta",
      shiny::span(class = "sb-badge", paste0("v", .sigbenchr_version)),
      shiny::span(
        class = "sb-badge sb-badge-warn",
        "Research use only — not a clinical decision tool"
      )
    )
  ),
  shiny::uiOutput("stale_note"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 3,
      shiny::div(
        class = "sb-actions",
        shiny::actionButton(
          "run", "Run analysis",
          class = "btn-primary btn-block",
          width = "100%"
        ),
        shiny::downloadButton(
          "download_report", "Download HTML report",
          class = "btn-default btn-block"
        ),
        shiny::downloadButton(
          "download_manifest", "Download manifest (JSON)",
          class = "btn-default btn-block"
        ),
        shiny::downloadButton(
          "download_pub_bundle", "Download publication bundle (.zip)",
          class = "btn-default btn-block"
        )
      ),
      shiny::tags$details(
        class = "sb-acc",
        open = NA,
        shiny::tags$summary("1. Signatures and cohorts"),
        shiny::checkboxGroupInput(
          "sig_choice", "Built-in signatures",
          choices = names(.builtin_sigs_ui),
          selected = intersect(
            c(
              "Cytokine-5", "IFNG6", "T_Cell_Inflamed_GEP",
              "Cytolytic_Activity", "Myeloid_Suppression"
            ),
            names(.builtin_sigs_ui)
          )
        ),
        shiny::textAreaInput(
          "custom_genes",
          "Custom signature (optional)",
          rows = 3,
          placeholder = "Separate genes by comma, space, or newline"
        ),
        shiny::tags$details(
          class = "sb-acc",
          shiny::tags$summary("Optional signed custom (genes_up / genes_down)"),
          shiny::textAreaInput(
            "custom_genes_up",
            "Up-regulated genes (genes_up)",
            rows = 2,
            placeholder = "e.g. IFNG, CXCL9, STAT1"
          ),
          shiny::textAreaInput(
            "custom_genes_down",
            "Down-regulated genes (genes_down)",
            rows = 2,
            placeholder = "e.g. VEGFA, TGFB1"
          )
        ),
        shiny::tags$p(
          class = "sb-hint",
          "Custom signatures appear as \"Custom\" in Jaccard and DeLong selectors; derived_from is omitted (independent = NA)."
        ),
        shiny::tags$p(
          class = "sb-hint",
          "Note: Bundled ICI cohorts contain ~2,800 genes (~14-15% coverage for external signatures). Load local CSV/TSV for full gene coverage."
        ),
        shiny::checkboxGroupInput(
          "datasets", "Cohorts (ICI)",
          choices = c(
            "van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019"
          ),
          selected = c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
        ),
        if (isTRUE(.sigbenchr_allow_upload)) {
          shiny::tags$details(
            class = "sb-acc",
            shiny::tags$summary("Upload local cohort (CSV/TSV)"),
            shiny::fileInput("upload_expr", "Expression matrix (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
            shiny::fileInput("upload_meta", "Clinical metadata (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
            shiny::checkboxInput("upload_genes_in_rows", "Genes in rows", TRUE),
            shiny::selectInput(
              "upload_gene_id_type", "Gene ID type",
              choices = c("HGNC", "Ensembl", "other"), selected = "HGNC"
            ),
            shiny::selectInput(
              "upload_normalization", "Normalization",
              choices = c("TPM", "FPKM", "normalized_counts", "log_normalized", "other"),
              selected = "TPM"
            ),
            shiny::selectInput(
              "upload_scale", "Scale",
              choices = c("linear", "log2"), selected = "log2"
            ),
            shiny::selectInput(
              "upload_species", "Species",
              choices = c("human", "mouse"), selected = "human"
            ),
            shiny::textInput("upload_sample_id_col", "sample_id column", value = "sample_id"),
            shiny::textInput("upload_response_col", "response column", value = "response_label"),
            shiny::textInput("upload_dataset_name", "Cohort name", value = "user_cohort"),
            shiny::actionButton("upload_add", "Validate and add cohort", class = "btn-info"),
            shiny::verbatimTextOutput("upload_status")
          )
        } else {
          shiny::tags$p(
            class = "sb-note",
            "Local cohort upload is enabled only via SigBenchR::run_app(allow_upload = TRUE)."
          )
        },
        shiny::tags$details(
          class = "sb-acc",
          shiny::tags$summary("Cohort registry"),
          shiny::selectInput(
            "registry_cohort", "Registered cohorts",
            choices = tryCatch(
              SigBenchR::list_cohort_registry()$cohort_id,
              error = function(e) c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
            )
          ),
          shiny::actionButton("registry_load", "Load from registry", class = "btn-default btn-sm"),
          shiny::verbatimTextOutput("registry_status")
        )
      ),
      shiny::tags$details(
        class = "sb-acc",
        open = NA,
        shiny::tags$summary("2. Method and clinical policy"),
        shiny::checkboxGroupInput(
          "methods", "Scoring methods",
          choices = c("mean", "zscore", "median", "rank", "ss_rank"),
          selected = c("mean", "zscore")
        ),
        shiny::radioButtons(
          "direction", "Signature direction",
          choices = c("High = responder (<)" = "<", "Low = responder (>)" = ">"),
          selected = "<"
        ),
        shiny::radioButtons(
          "sd_handling", "SD policy",
          choices = c(
            "Exclude SD [CR/PR vs PD]" = "exclude",
            "ORR [CR/PR vs SD/PD]" = "as_0",
            "DCR [CR/PR/SD vs PD]" = "as_1"
          ),
          selected = "exclude"
        ),
        shiny::checkboxInput(
          "compute_null",
          "Null permutation (recommended random gene-set control)",
          FALSE
        ),
        shiny::tags$p(
          class = "sb-hint",
          "Recommended: compare signature AUC against same-sized random gene sets (compute_null)."
        )
      ),
      shiny::tags$details(
        class = "sb-acc",
        shiny::tags$summary("3. Advanced statistical settings"),
        shiny::sliderInput("coverage_threshold", "Coverage threshold", 0.5, 1, 0.70, 0.05),
        shiny::sliderInput("sensitivity_threshold", "Robustness threshold (arbitrary)", 0.05, 0.3, 0.10, 0.01),
        shiny::checkboxInput("compute_pr", "PR-AUC + bootstrap CI", TRUE),
        shiny::numericInput("n_boot", "PR bootstrap (n_boot)", value = 200, min = 0, max = 2000, step = 50),
        shiny::numericInput("n_perm", "Null B", value = 1000, min = 10, max = 2000, step = 50),
        shiny::checkboxInput(
          "match_expression",
          "Null: expression-matched sampling",
          FALSE
        ),
        shiny::numericInput("seed", "Seed", value = 42, min = 1, step = 1)
      )
    ),
    shiny::mainPanel(
      width = 9,
      shiny::uiOutput("kpi_strip"),
      shiny::tabsetPanel(
        id = "main_tabs",
        shiny::tabPanel(
          "Overview",
          shiny::div(
            class = "sb-panel-block",
            shiny::radioButtons(
              "heatmap_method", "Active view method",
              choices = c("mean", "zscore"), selected = "mean",
              inline = TRUE
            ),
            shiny::tags$h4(class = "sb-section", "AUC heatmap"),
            shiny::plotOutput("heatmap", height = "320px"),
            shiny::uiOutput("overview_stat_alerts"),
            shiny::tags$h4(class = "sb-section", "Results table (ROC-AUC, PR-AUC & Null Permutation FDR)"),
            shiny::tags$details(
              class = "sb-acc",
              shiny::tags$summary("Statistical column glossary (raw vs BH-adjusted p-values, PR diagnostics & background universe)"),
              shiny::tags$ul(
                shiny::tags$li(shiny::tags$strong("auc / auc_ci_low / auc_ci_high:"), " DeLong ROC-AUC point estimate and 95% confidence interval."),
                shiny::tags$li(shiny::tags$strong("pr_auc / pr_ci_low / pr_ci_high / pr_baseline:"), " Precision-Recall AUC, stratified bootstrap 95% CI, and positive-class prevalence baseline (n_pos / n)."),
                shiny::tags$li(shiny::tags$strong("small_positive_class / ci_unreliable / n_boot_valid:"), " PR-AUC diagnostic flags: warns when n_pos < 5 (or n_neg < 5) or when tie fraction / valid bootstrap replicates make the PR bootstrap CI unstable."),
                shiny::tags$li(shiny::tags$strong("null_p (Raw) vs null_p_adj_bh (Within-Method BH) vs null_p_adj_bh_global (Global BH):"), " Empirical random-gene permutation p-value (raw), Benjamini-Hochberg FDR adjusted within each scoring method, and BH-adjusted globally across all tested methods x signatures x cohorts."),
                shiny::tags$li(shiny::tags$strong("background_type:"), " Null sampling universe classification (Subsampled panel cohort (~2,800 genes), Full transcriptome (>=5,000 genes), or Partial ortholog-mapped subset).")
              )
            ),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("table")
            } else {
              shiny::tableOutput("table_fallback")
            }
          )
        ),
        shiny::tabPanel(
          "ROC and Forest",
          shiny::div(
            class = "sb-panel-block",
            shiny::div(
              class = "sb-inline-controls",
              shiny::selectInput(
                "roc_dataset", "Cohort",
                choices = c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
              )
            ),
            shiny::tags$h4(class = "sb-section", "ROC curves"),
            shiny::plotOutput("roc_plot", height = "420px")
          ),
          shiny::div(
            class = "sb-panel-block",
            shiny::div(
              class = "sb-inline-controls",
              shiny::selectInput(
                "forest_sig", "Forest signature",
                choices = c("Cytokine-5", "IFNG6", "Random3")
              )
            ),
            shiny::tags$h4(class = "sb-section", "AUC forest"),
            shiny::plotOutput("forest_plot", height = "320px")
          )
        ),
        shiny::tabPanel(
          "Summary and subgroups",
          shiny::div(
            class = "sb-panel-block",
            shiny::checkboxInput(
              "stratify_by",
              "Stratify by treatment / endpoint (stratify_by)",
              TRUE
            ),
            shiny::checkboxInput(
              "include_nonindependent",
              "Include Custom / non-independent signatures in summary",
              FALSE
            ),
            shiny::tags$p(
              class = "sb-hint",
              "When checked, summary uses stratify_by = c(\"treatment\", \"endpoint_type\") (e.g. anti-CTLA4 clinical benefit vs anti-PD1 RECIST)."
            ),
            shiny::tags$h4(class = "sb-section", "Weighted summary (n + inverse-variance, Q, I2)"),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("summary_table")
            } else {
              shiny::tableOutput("summary_fallback")
            }
          ),
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "Scoring robustness (check_sensitivity)"),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("sensitivity_table")
            } else {
              shiny::tableOutput("sensitivity_fallback")
            }
          ),
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "Clinical survival benchmark (Cox HR & Log-rank FDR)"),
            shiny::tags$p(
              class = "sb-hint",
              "Evaluated on cohorts containing time-to-event columns (os_time / os_event). Continuous Cox PH (p_adj_bh) and median-dichotomized Log-rank (logrank_p_adj_bh) are BH-adjusted as separate hypothesis families; Schoenfeld PH assumption test (ph_p) is left unadjusted as a model diagnostic."
            ),
            shiny::uiOutput("survival_alerts"),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("survival_table")
            } else {
              shiny::tableOutput("survival_fallback")
            }
          )
        ),
        shiny::tabPanel(
          "Comparison",
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "Paired DeLong"),
            shiny::div(
              class = "sb-inline-controls",
              shiny::selectInput(
                "cmp_dataset", "Dataset",
                choices = c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
              ),
              shiny::selectInput(
                "cmp_sig1", "Signature 1",
                choices = c("Cytokine-5", "IFNG6", "Random3")
              ),
              shiny::selectInput(
                "cmp_sig2", "Signature 2",
                choices = c("IFNG6", "Cytokine-5", "Random3")
              ),
              shiny::actionButton("run_cmp", "Compare", class = "btn-primary")
            ),
            shiny::verbatimTextOutput("cmp_status"),
            shiny::tags$p(
              class = "sb-hint",
              "For a single comparison, p_adj_bh equals the raw p-value (p_adj_bh == p_value)."
            ),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("cmp_table")
            } else {
              shiny::tableOutput("cmp_fallback")
            }
          ),
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "Jaccard overlap (including Custom)"),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("jaccard_table")
            } else {
              shiny::tableOutput("jaccard_fallback")
            }
          ),
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "AUC correlation (across datasets)"),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("corr_table")
            } else {
              shiny::tableOutput("corr_fallback")
            }
          )
        ),
        shiny::tabPanel(
          "Gene Contributions (LOO)",
          shiny::div(
            class = "sb-panel-block",
            shiny::div(
              class = "sb-inline-controls",
              shiny::selectInput(
                "gene_contrib_sig", "Signature",
                choices = c("Cytokine-5", "IFNG6", "Random3")
              )
            ),
            shiny::tags$p(
              class = "sb-hint",
              "Leave-one-out delta-AUC = full_auc - loo_auc across cohorts in the last run snapshot. Positive delta-AUC identifies signal drivers."
            ),
            shiny::plotOutput("gene_contrib_plot", height = "340px"),
            shiny::tags$h4(class = "sb-section", "Per-gene LOO summary & immune axes"),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("gene_contrib_table")
            } else {
              shiny::tableOutput("gene_contrib_fallback")
            }
          )
        ),
        shiny::tabPanel(
          "Beyond Infiltration",
          shiny::div(
            class = "sb-panel-block",
            shiny::div(
              class = "sb-inline-controls",
              shiny::selectInput(
                "beyond_inf_sig", "Signature",
                choices = c("Cytokine-5", "IFNG6", "Random3")
              ),
              shiny::selectInput(
                "beyond_inf_axis", "Reference infiltration axis",
                choices = c("CD8_T_Cell", "IFN_Gamma", "Cytotoxicity"),
                selected = "CD8_T_Cell"
              )
            ),
            shiny::tags$p(
              class = "sb-hint",
              "Tests whether the signature adds predictive value beyond baseline T-cell infiltration via likelihood-ratio test (glm: response ~ sig_score_z + inf_score_z vs inf_score_z) and fold-inside paired cross-validated delta-AUC."
            ),
            shiny::uiOutput("beyond_inf_alerts"),
            shiny::verbatimTextOutput("beyond_inf_summary"),
            shiny::tags$h4(class = "sb-section", "Fold-level paired CV delta-AUC (+/- 95% CI)"),
            shiny::plotOutput("beyond_inf_plot", height = "280px"),
            shiny::tags$details(
              class = "sb-acc",
              shiny::tags$summary("Diagnostic field explanations (EPV, VIF, BH-adjusted LRT, Stouffer Z & analysis status)"),
              shiny::tags$ul(
                shiny::tags$li(shiny::tags$strong("lrt_p (Raw) vs lrt_p_adj_bh (BH-adjusted FDR):"), " Likelihood-ratio test p-value before and after Benjamini-Hochberg correction across evaluated cohorts. Independent value requires lrt_p_adj_bh < 0.05."),
                shiny::tags$li(shiny::tags$strong("delta_auc_combined / sd_delta_auc / [delta_auc_ci_low, delta_auc_ci_high]:"), " Paired cross-validated test AUC gain of the combined model over infiltration alone, with fold-level standard deviation and 95% confidence interval (preprocessing fit strictly inside each training fold)."),
                shiny::tags$li(shiny::tags$strong("epv (Events Per Variable) & low_epv:"), " min(n_pos, n_neg) / 2 for the 2-predictor logistic model. When EPV < 5, maximum-likelihood estimates are exploratory."),
                shiny::tags$li(shiny::tags$strong("vif (Variance Inflation Factor) & high_collinearity:"), " 1 / (1 - r^2) between signature and reference infiltration axis. VIF >= 5 warns of multicollinearity; VIF >= 10 aborts coefficient estimation as collinearity failure."),
                shiny::tags$li(shiny::tags$strong("direction_conflict & combined_lrt_p (Directional Stouffer Z):"), " Cross-cohort synthesis uses sample-size weighted directional Stouffer Z so opposite-sign effects cancel rather than inflate significance.")
              )
            ),
            shiny::tags$h4(class = "sb-section", "Per-cohort adjusted logistic regression"),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("beyond_inf_table")
            } else {
              shiny::tableOutput("beyond_inf_fallback")
            }
          )
        ),
        shiny::tabPanel(
          "Coverage and reproducibility",
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "Gene coverage and missing genes"),
            shiny::tags$p(
              class = "sb-hint",
              "Per signature x cohort: genes found/total, coverage percent, and missing_genes."
            ),
            if (requireNamespace("DT", quietly = TRUE)) {
              DT::DTOutput("coverage_table")
            } else {
              shiny::tableOutput("coverage_fallback")
            }
          ),
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "Copyable R code"),
            shiny::tags$p(
              class = "sb-hint",
              "Generated from the last run snapshot parameters; does not drift with live inputs."
            ),
            shiny::textAreaInput(
              "r_code_text",
              label = NULL,
              value = "",
              rows = 16,
              width = "100%"
            ),
            shiny::tags$button(
              "Copy code",
              id = "copy_r_code",
              class = "btn btn-default",
              type = "button",
              onclick = "sigbenchrCopyRcode();"
            )
          ),
          shiny::div(
            class = "sb-panel-block",
            shiny::tags$h4(class = "sb-section", "Scientific limits and provenance"),
            shiny::tags$ul(
              shiny::tags$li("Research exploration tool only; not designed for clinical decisions."),
              shiny::tags$li("Uses public melanoma ICI cohorts (datasets_ici) and demo signatures."),
              shiny::tags$li("Plots, R code, and HTML report are driven only by the results() snapshot (res, sigs, params)."),
              shiny::tags$li("Cytokine-5 includes derived_from (independent = FALSE); IFNG6 / Random3 are independent; Custom is independent = NA.")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Shiny varsayilani gercek ICI kohortlari; sentetik demo_A/demo_B testlerde kalir.
  # Kullanici kohortlari (yalniz allow_upload = TRUE) reactiveVal ile eklenir;
  # snapshot kilidi current_params / results ile korunur.
  user_datasets <- shiny::reactiveVal(list())
  upload_status_msg <- shiny::reactiveVal("")
  registry_status_msg <- shiny::reactiveVal("")

  output$upload_status <- shiny::renderText(upload_status_msg())
  output$registry_status <- shiny::renderText(registry_status_msg())

  active_dataset_pool <- shiny::reactive({
    c(datasets_ici, user_datasets())
  })

  shiny::observeEvent(input$upload_add, {
    shiny::req(isTRUE(.sigbenchr_allow_upload))
    shiny::req(input$upload_expr, input$upload_meta)
    nm <- trimws(as.character(input$upload_dataset_name))
    if (!nzchar(nm)) nm <- "user_cohort"
    warn_msgs <- character(0)
    ds <- tryCatch(
      withCallingHandlers(
        SigBenchR::read_user_dataset(
          expr_path = input$upload_expr$datapath,
          meta_path = input$upload_meta$datapath,
          genes_in_rows = isTRUE(input$upload_genes_in_rows),
          gene_id_type = as.character(input$upload_gene_id_type),
          normalization = as.character(input$upload_normalization),
          scale = as.character(input$upload_scale),
          sample_id_col = as.character(input$upload_sample_id_col),
          response_col = as.character(input$upload_response_col),
          dataset_name = nm
        ),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    if (inherits(ds, "error")) {
      upload_status_msg(paste0("Error: ", conditionMessage(ds)))
      return()
    }
    ds$info$species <- if (!is.null(input$upload_species) && nzchar(input$upload_species)) {
      as.character(input$upload_species)
    } else {
      "human"
    }
    pool <- user_datasets()
    pool[[nm]] <- ds
    user_datasets(pool)
    all_names <- c(names(datasets_ici), names(pool))
    cur_sel <- as.character(input$datasets)
    cur_sel <- cur_sel[nzchar(cur_sel)]
    shiny::updateCheckboxGroupInput(
      session, "datasets",
      choices = all_names,
      selected = unique(c(cur_sel, nm))
    )
    msg <- paste0(
      "Added: ", nm,
      " | genes=", ds$info$n_genes,
      " | samples=", ds$info$n_samples,
      " | norm=", ds$info$normalization,
      " | scale=", ds$info$scale
    )
    if (length(warn_msgs) > 0L) {
      msg <- paste0(msg, "\nWarnings:\n- ", paste(warn_msgs, collapse = "\n- "))
    }
    upload_status_msg(msg)
  })

  shiny::observeEvent(input$registry_load, {
    shiny::req(input$registry_cohort)
    cid <- as.character(input$registry_cohort)[1]
    ds <- tryCatch(
      SigBenchR::load_registered_cohort(cid),
      error = function(e) e
    )
    if (inherits(ds, "error")) {
      registry_status_msg(paste0("Error: ", conditionMessage(ds)))
      return()
    }
    pool <- user_datasets()
    if (!cid %in% names(datasets_ici)) {
      pool[[cid]] <- ds
      user_datasets(pool)
    }
    all_names <- unique(c(names(datasets_ici), names(pool)))
    cur_sel <- as.character(input$datasets)
    cur_sel <- cur_sel[nzchar(cur_sel)]
    shiny::updateCheckboxGroupInput(
      session, "datasets",
      choices = all_names,
      selected = unique(c(cur_sel, cid))
    )
    registry_status_msg(sprintf(
      "Loaded '%s' (species=%s, genes=%d, samples=%d)",
      cid, ds$info$species, nrow(ds$expr), ncol(ds$expr)
    ))
  })

  active_sig_names <- shiny::reactive({
    chosen <- as.character(input$sig_choice)
    chosen <- chosen[nzchar(chosen)]
    custom_txt <- if (is.null(input$custom_genes)) "" else as.character(input$custom_genes)[1]
    up_txt <- if (is.null(input$custom_genes_up)) "" else as.character(input$custom_genes_up)[1]
    down_txt <- if (is.null(input$custom_genes_down)) "" else as.character(input$custom_genes_down)[1]
    if (nzchar(trimws(custom_txt)) || nzchar(trimws(up_txt)) || nzchar(trimws(down_txt))) {
      chosen <- c(chosen, "Custom")
    }
    unique(chosen)
  })

  build_sigs <- function(direction,
                         sig_choice,
                         custom_genes,
                         custom_genes_up = NULL,
                         custom_genes_down = NULL) {
    chosen <- as.character(sig_choice)
    chosen <- chosen[nzchar(chosen)]
    sigs <- list()
    for (nm in chosen) {
      shiny::validate(shiny::need(nm %in% names(.builtin_sigs_pool), paste0("Unknown signature: ", nm)))
      sig <- .builtin_sigs_pool[[nm]]
      if (is.null(sig$direction) || !nzchar(as.character(sig$direction)[1])) {
        sig$direction <- direction
      }
      sigs[[nm]] <- sig
    }
    custom_txt <- if (is.null(custom_genes)) "" else as.character(custom_genes)[1]
    up_txt <- if (is.null(custom_genes_up)) "" else as.character(custom_genes_up)[1]
    down_txt <- if (is.null(custom_genes_down)) "" else as.character(custom_genes_down)[1]
    has_flat <- !is.na(custom_txt) && nzchar(trimws(custom_txt))
    has_updown <- (!is.na(up_txt) && nzchar(trimws(up_txt))) ||
      (!is.na(down_txt) && nzchar(trimws(down_txt)))

    shiny::validate(shiny::need(
      !(has_flat && has_updown),
      "Provide either 'Custom signature' or signed 'genes_up / genes_down', not both."
    ))

    if (has_updown) {
      v_res <- tryCatch(
        .validate_custom_signature(
          genes_up = up_txt,
          genes_down = down_txt,
          custom_name = "Custom",
          existing_names = names(sigs)
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(
        !inherits(v_res, "error"),
        if (inherits(v_res, "error")) conditionMessage(v_res) else ""
      ))
      # derived_from omitted -> independent = NA (fail-closed)
      sigs[[v_res$name]] <- list(
        genes_up = if (length(v_res$genes_up) > 0L) v_res$genes_up else NULL,
        genes_down = if (length(v_res$genes_down) > 0L) v_res$genes_down else NULL,
        direction = direction
      )
    } else if (has_flat) {
      v_res <- tryCatch(
        .validate_custom_signature(
          genes_text = custom_txt,
          custom_name = "Custom",
          existing_names = names(sigs)
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(
        !inherits(v_res, "error"),
        if (inherits(v_res, "error")) conditionMessage(v_res) else ""
      ))
      # derived_from omitted -> independent = NA (fail-closed)
      sigs[[v_res$name]] <- list(genes = v_res$genes, direction = direction)
    }
    shiny::validate(shiny::need(length(sigs) >= 1, "Select at least one built-in signature or add a custom signature."))
    sigs
  }

  .safe_int <- function(val) {
    if (is.null(val) || length(val) == 0L) return(NA_integer_)
    v <- suppressWarnings(as.integer(val[1]))
    if (length(v) == 0L) NA_integer_ else v
  }

  current_params <- shiny::reactive({
    hm_cur <- as.character(input$heatmap_method)
    hm_cur <- hm_cur[nzchar(hm_cur)]
    if (length(hm_cur) < 1L) {
      m_cur <- as.character(input$methods)
      m_cur <- m_cur[nzchar(m_cur)]
      hm_cur <- if (length(m_cur) >= 1L) m_cur[1] else character(0)
    }
    fs_cur <- as.character(input$forest_sig)
    fs_cur <- fs_cur[nzchar(fs_cur)]
    if (length(fs_cur) < 1L) {
      s_cur <- active_sig_names()
      fs_cur <- if (length(s_cur) >= 1L) s_cur[1] else character(0)
    }
    list(
      methods = as.character(input$methods),
      direction = as.character(input$direction),
      coverage_threshold = as.numeric(input$coverage_threshold),
      sensitivity_threshold = as.numeric(input$sensitivity_threshold),
      n_boot = .safe_int(input$n_boot),
      n_perm = .safe_int(input$n_perm),
      seed = .safe_int(input$seed),
      compute_pr = isTRUE(input$compute_pr),
      compute_null = isTRUE(input$compute_null),
      sd_handling = as.character(input$sd_handling),
      match_expression = isTRUE(input$match_expression),
      datasets = as.character(input$datasets),
      sig_choice = as.character(input$sig_choice),
      custom_genes = if (is.null(input$custom_genes)) "" else as.character(input$custom_genes),
      custom_genes_up = if (is.null(input$custom_genes_up)) "" else as.character(input$custom_genes_up),
      custom_genes_down = if (is.null(input$custom_genes_down)) "" else as.character(input$custom_genes_down),
      heatmap_method = hm_cur,
      forest_signature = fs_cur
    )
  })

  # Dinamik dropdown: yontemler
  shiny::observe({
    meth <- as.character(input$methods)
    meth <- meth[nzchar(meth)]
    if (length(meth) < 1L) {
      return()
    }
    cur <- input$heatmap_method
    sel <- if (!is.null(cur) && cur %in% meth) cur else meth[1]
    shiny::updateRadioButtons(
      session, "heatmap_method",
      choices = meth,
      selected = sel
    )
  })

  # Dinamik dropdown: imzalar (Custom dahil)
  shiny::observe({
    sigs <- active_sig_names()
    if (length(sigs) < 1L) {
      sigs <- c("Cytokine-5", "IFNG6")
    }
    cur_f <- input$forest_sig
    sel_f <- if (!is.null(cur_f) && cur_f %in% sigs) cur_f else sigs[1]
    shiny::updateSelectInput(session, "forest_sig", choices = sigs, selected = sel_f)

    cur1 <- input$cmp_sig1
    cur2 <- input$cmp_sig2
    sel1 <- if (!is.null(cur1) && cur1 %in% sigs) cur1 else sigs[1]
    sel2 <- if (!is.null(cur2) && cur2 %in% sigs) {
      cur2
    } else if (length(sigs) >= 2L) {
      sigs[2]
    } else {
      sigs[1]
    }
    shiny::updateSelectInput(session, "cmp_sig1", choices = sigs, selected = sel1)
    shiny::updateSelectInput(session, "cmp_sig2", choices = sigs, selected = sel2)
    cur_gc <- input$gene_contrib_sig
    sel_gc <- if (!is.null(cur_gc) && cur_gc %in% sigs) cur_gc else sigs[1]
    shiny::updateSelectInput(session, "gene_contrib_sig", choices = sigs, selected = sel_gc)
    cur_bi <- input$beyond_inf_sig
    sel_bi <- if (!is.null(cur_bi) && cur_bi %in% sigs) cur_bi else sigs[1]
    shiny::updateSelectInput(session, "beyond_inf_sig", choices = sigs, selected = sel_bi)
  })

  # Secili kohortlara gore ROC / DeLong dataset secenekleri
  shiny::observeEvent(input$datasets, {
    ds <- as.character(input$datasets)
    ds <- ds[nzchar(ds)]
    if (length(ds) < 1L) {
      return()
    }
    cur_roc <- input$roc_dataset
    cur_cmp <- input$cmp_dataset
    sel_roc <- if (!is.null(cur_roc) && cur_roc %in% ds) cur_roc else ds[1]
    sel_cmp <- if (!is.null(cur_cmp) && cur_cmp %in% ds) cur_cmp else ds[1]
    shiny::updateSelectInput(session, "roc_dataset", choices = ds, selected = sel_roc)
    shiny::updateSelectInput(session, "cmp_dataset", choices = ds, selected = sel_cmp)
  }, ignoreNULL = FALSE)

  results <- shiny::eventReactive(input$run, {
    shiny::req(length(input$datasets) > 0, length(input$methods) > 0)
    params <- current_params()
    min_boot <- if (isTRUE(params$compute_pr)) 1L else 0L
    shiny::validate(
      shiny::need(
        !is.na(params$n_boot) && is.finite(params$n_boot) && params$n_boot >= min_boot,
        paste0("n_boot must be an integer >= ", min_boot, ".")
      ),
      shiny::need(
        !is.na(params$n_perm) && is.finite(params$n_perm) && params$n_perm >= 1L,
        "n_perm must be an integer >= 1."
      ),
      shiny::need(
        !is.na(params$seed) && is.finite(params$seed) && params$seed >= 1L,
        "Seed must be an integer >= 1."
      )
    )
    sigs <- build_sigs(
      params$direction,
      params$sig_choice,
      params$custom_genes,
      params$custom_genes_up,
      params$custom_genes_down
    )
    pool <- active_dataset_pool()
    miss_ds <- setdiff(params$datasets, names(pool))
    shiny::validate(shiny::need(
      length(miss_ds) == 0L,
      paste0("Unknown dataset: ", paste(miss_ds, collapse = ", "))
    ))
    ds <- pool[params$datasets]
    hm_init <- as.character(params$heatmap_method)
    hm_init <- hm_init[nzchar(hm_init)]
    params$heatmap_method <- if (length(hm_init) >= 1L && hm_init[1] %in% params$methods) {
      hm_init[1]
    } else {
      params$methods[1]
    }
    fs_init <- as.character(params$forest_signature)
    fs_init <- fs_init[nzchar(fs_init)]
    params$forest_signature <- if (length(fs_init) >= 1L && fs_init[1] %in% names(sigs)) {
      fs_init[1]
    } else {
      names(sigs)[1]
    }

    shiny::withProgress(message = "Running benchmark...", value = 0.2, {
      shiny::incProgress(0.3, detail = "ROC / PR / null")
      out <- SigBenchR::run_benchmark(
        signatures = sigs,
        datasets = ds,
        methods = params$methods,
        direction = params$direction,
        coverage_threshold = params$coverage_threshold,
        sd_handling = params$sd_handling,
        compute_pr = params$compute_pr,
        n_boot = params$n_boot,
        compute_null = params$compute_null,
        n_perm = params$n_perm,
        seed = params$seed,
        match_expression = params$match_expression
      )
      shiny::incProgress(0.5)
      list(
        res = out,
        sigs = sigs,
        datasets = ds,
        params = params
      )
    })
  }, ignoreInit = TRUE)

  stale_inputs <- shiny::reactive({
    if (is.null(input$run) || identical(as.integer(input$run), 0L)) {
      return(FALSE)
    }
    snap <- results()
    if (is.null(snap)) {
      return(FALSE)
    }
    cur <- current_params()
    keys <- c(
      "methods", "direction", "coverage_threshold", "sensitivity_threshold",
      "n_boot", "n_perm", "seed",
      "compute_pr", "compute_null", "sd_handling", "match_expression",
      "datasets", "sig_choice", "custom_genes",
      "custom_genes_up", "custom_genes_down",
      "heatmap_method", "forest_signature"
    )
    !identical(cur[keys], snap$params[keys])
  })

  output$stale_note <- shiny::renderUI({
    if (isTRUE(stale_inputs())) {
      shiny::div(
        class = "sb-stale",
        "Parameters changed — press 'Run analysis' to refresh results."
      )
    } else {
      NULL
    }
  })

  report_model <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    hm_snap <- .snap_method(snap)
    fsig_snap <- .snap_forest_sig(snap)
    snap_params <- utils::modifyList(
      snap$params,
      list(heatmap_method = hm_snap, forest_signature = fsig_snap)
    )
    bio_ctx <- tryCatch(
      SigBenchR::correlate_reference_axes(
        datasets = snap$datasets,
        signatures = snap$sigs,
        method = hm_snap,
        sd_handling = snap$params$sd_handling
      ),
      error = function(e) NULL
    )
    surv_ctx <- tryCatch(
      suppressWarnings(
        SigBenchR::run_survival_benchmark(
          datasets = snap$datasets,
          signatures = snap$sigs,
          method = hm_snap
        )
      ),
      error = function(e) NULL
    )
    bi_ctx <- tryCatch(
      suppressWarnings({
        bi_list <- list()
        for (sn in names(snap$sigs)) {
          bi_list[[sn]] <- SigBenchR::test_beyond_infiltration(
            datasets = snap$datasets,
            signature = snap$sigs[[sn]],
            infiltration_axis = "CD8_T_Cell",
            method = hm_snap,
            sd_handling = snap$params$sd_handling,
            direction = snap$params$direction
          )
        }
        bi_list
      }),
      error = function(e) NULL
    )
    SigBenchR::build_report_model(
      results = snap$res,
      summary_tbl = NULL,
      sensitivity_tbl = sens(),
      params_list = snap_params,
      ref_axes = bio_ctx,
      surv_res = surv_ctx,
      beyond_inf = bi_ctx
    )
  })

  output$kpi_strip <- shiny::renderUI({
    if (is.null(input$run) || identical(as.integer(input$run), 0L)) {
      return(shiny::div(
        class = "sb-welcome",
        shiny::tags$h4("Getting started"),
        shiny::tags$ol(
          shiny::tags$li("Select signatures and cohorts on the left (add a custom gene list if needed)."),
          shiny::tags$li("Set scoring method, signature direction, and SD policy."),
          shiny::tags$li("Press \"Run analysis\"; KPI cards and tabs fill from the snapshot."),
          shiny::tags$li("Inspect ROC / Forest, summary subgroups, and DeLong comparisons in their tabs.")
        )
      ))
    }
    snap <- results()
    if (is.null(snap)) {
      return(shiny::div(
        class = "sb-welcome",
        shiny::tags$h4("Getting started"),
        shiny::tags$ol(
          shiny::tags$li("Select signatures and cohorts on the left (add a custom gene list if needed)."),
          shiny::tags$li("Set scoring method, signature direction, and SD policy."),
          shiny::tags$li("Press \"Run analysis\"; KPI cards and tabs fill from the snapshot."),
          shiny::tags$li("Inspect ROC / Forest, summary subgroups, and DeLong comparisons in their tabs.")
        )
      ))
    }

    m <- report_model()
    best <- m$tier1$best_signature
    alerts <- m$tier1$alerts

    best_label <- if (!is.null(best) && nrow(best) > 0L) {
      sprintf("%s (%s, k=%d)", best$signature[1], best$method[1], best$k[1])
    } else {
      "No independent result"
    }
    pooled_auc_label <- if (!is.null(best) && nrow(best) > 0L && is.finite(best$auc_re[1])) {
      if (is.finite(best$auc_re_ci_low[1]) && is.finite(best$auc_re_ci_high[1])) {
        sprintf("%.3f (%.3f\u2013%.3f)", best$auc_re[1], best$auc_re_ci_low[1], best$auc_re_ci_high[1])
      } else {
        sprintf("%.3f", best$auc_re[1])
      }
    } else {
      "\u2014"
    }
    i2_label <- if (!is.null(best) && nrow(best) > 0L && is.finite(best$i2_logit[1])) {
      sprintf("%.1f%%", best$i2_logit[1])
    } else {
      "\u2014"
    }
    tier_txt <- if (!is.null(best) && nrow(best) > 0L) best$priority_tier[1] else "Unevaluable"
    tier_cls <- if (grepl("high", tolower(tier_txt))) {
      "sb-tier-pill sb-tier-high"
    } else if (grepl("medium", tolower(tier_txt))) {
      "sb-tier-pill sb-tier-med"
    } else {
      "sb-tier-pill sb-tier-low"
    }

    alert_ui <- if (is.data.frame(alerts) && nrow(alerts) > 0L) {
      shiny::div(
        class = "sb-alert-strip",
        lapply(seq_len(nrow(alerts)), function(i) {
          cls <- if (identical(alerts$severity[i], "info")) "sb-alert-item sb-alert-info" else "sb-alert-item"
          shiny::div(
            class = cls,
            shiny::tags$strong(paste0("[", alerts$rule_id[i], "] ", alerts$title[i], ": ")),
            alerts$message[i]
          )
        })
      )
    } else {
      NULL
    }

    n_eval <- sum(!(snap$res$skipped %in% TRUE), na.rm = TRUE)
    n_skip <- sum(snap$res$skipped %in% TRUE, na.rm = TRUE)
    stage_bar <- shiny::div(
      class = "sb-stage-bar",
      shiny::span(
        shiny::tags$strong("Active snapshot: "),
        sprintf(
          "%d signature(s) x %d cohort(s) x %d method(s) | Evaluated: %d | Skipped: %d",
          length(snap$sigs), length(snap$datasets), length(snap$params$methods), n_eval, n_skip
        )
      ),
      shiny::span(
        sprintf(
          "PR-AUC: %s | Null permutation: %s (seed=%d)",
          if (isTRUE(snap$params$compute_pr)) paste0("ON (B=", snap$params$n_boot, ")") else "OFF",
          if (isTRUE(snap$params$compute_null)) paste0("ON (B=", snap$params$n_perm, ")") else "OFF",
          snap$params$seed
        )
      )
    )

    shiny::tagList(
      stage_bar,
      shiny::div(
        class = "sb-kpi-row",
        shiny::div(
          class = "sb-kpi",
          shiny::tags$p(class = "sb-kpi-label", "Best independent signature"),
          shiny::tags$p(class = "sb-kpi-value", best_label)
        ),
        shiny::div(
          class = "sb-kpi",
          shiny::tags$p(class = "sb-kpi-label", "Pooled REML AUC (95% CI)"),
          shiny::tags$p(class = "sb-kpi-value", pooled_auc_label)
        ),
        shiny::div(
          class = "sb-kpi",
          shiny::tags$p(class = "sb-kpi-label", "Logit heterogeneity (I2 logit)"),
          shiny::tags$p(class = "sb-kpi-value", i2_label)
        ),
        shiny::div(
          class = "sb-kpi",
          shiny::tags$p(class = "sb-kpi-label", "Evidence tier"),
          shiny::tags$p(class = "sb-kpi-value", shiny::span(class = tier_cls, tier_txt))
        )
      ),
      alert_ui
    )
  })

  sens <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    SigBenchR::check_sensitivity(
      snap$res,
      threshold = snap$params$sensitivity_threshold
    )
  })

  .validated_hm_method <- function(snap) {
    hm <- if (!is.null(input$heatmap_method) &&
              length(input$heatmap_method) >= 1L &&
              nzchar(as.character(input$heatmap_method)[1])) {
      as.character(input$heatmap_method)[1]
    } else {
      as.character(snap$params$heatmap_method)[1]
    }
    shiny::validate(shiny::need(
      !is.na(hm) && nzchar(hm) && hm %in% snap$params$methods,
      paste0("Selected method '", hm, "' is not in the last run snapshot; press 'Run analysis'.")
    ))
    hm
  }

  .snap_method <- function(snap) {
    hm <- as.character(snap$params$heatmap_method)
    hm <- hm[nzchar(hm)]
    if (length(hm) >= 1L && hm[1] %in% snap$params$methods) {
      hm[1]
    } else {
      snap$params$methods[1]
    }
  }

  .snap_forest_sig <- function(snap) {
    fs <- as.character(snap$params$forest_signature)
    fs <- fs[nzchar(fs)]
    if (length(fs) >= 1L && fs[1] %in% names(snap$sigs)) {
      fs[1]
    } else {
      names(snap$sigs)[1]
    }
  }

  summary_res <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    hm <- .validated_hm_method(snap)
    strat <- if (isTRUE(input$stratify_by)) {
      c("treatment", "endpoint_type")
    } else {
      NULL
    }
    suppressWarnings(
      SigBenchR::summarize_auc(
        snap$res,
        method = hm,
        include_nonindependent = isTRUE(input$include_nonindependent),
        stratify_by = strat
      )
    )
  })

  # Unstratified summary for forest reference (avoids multi-row stratify_by clash)
  summary_res_forest <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    hm <- .validated_hm_method(snap)
    suppressWarnings(
      SigBenchR::summarize_auc(snap$res, method = hm, stratify_by = NULL)
    )
  })

  shiny::observe({
    if (is.null(input$run) || identical(as.integer(input$run), 0L)) {
      shiny::updateTextAreaInput(
        session, "r_code_text",
        value = "# Run the analysis to generate reproducible R code here."
      )
      return()
    }
    snap <- results()
    if (is.null(snap)) {
      shiny::updateTextAreaInput(
        session, "r_code_text",
        value = "# Run the analysis to generate reproducible R code here."
      )
      return()
    }
    p <- snap$params
    demo_names <- as.character(p$sig_choice)
    demo_names <- demo_names[nzchar(demo_names)]
    custom_txt <- if (is.null(p$custom_genes)) "" else as.character(p$custom_genes)[1]
    up_txt <- if (is.null(p$custom_genes_up)) "" else as.character(p$custom_genes_up)[1]
    down_txt <- if (is.null(p$custom_genes_down)) "" else as.character(p$custom_genes_down)[1]
    genes <- if (nzchar(trimws(custom_txt))) {
      SigBenchR::parse_genes(custom_txt)
    } else {
      NULL
    }
    genes_up <- if (nzchar(trimws(up_txt))) {
      SigBenchR::parse_genes(up_txt)
    } else {
      NULL
    }
    genes_down <- if (nzchar(trimws(down_txt))) {
      SigBenchR::parse_genes(down_txt)
    } else {
      NULL
    }
    sig_arg <- if (length(demo_names) > 0L) demo_names else NULL
    code <- tryCatch(
      SigBenchR::generate_r_code(
        signature_name = sig_arg,
        genes = genes,
        methods = p$methods,
        direction = p$direction,
        coverage_threshold = p$coverage_threshold,
        n_boot = p$n_boot,
        n_perm = p$n_perm,
        seed = p$seed,
        dataset_names = p$datasets,
        compute_pr = p$compute_pr,
        compute_null = p$compute_null,
        sd_handling = p$sd_handling,
        match_expression = p$match_expression,
        genes_up = genes_up,
        genes_down = genes_down
      ),
      error = function(e) paste0("# Could not generate code: ", conditionMessage(e))
    )
    shiny::updateTextAreaInput(session, "r_code_text", value = code)
  })

  output$download_report <- shiny::downloadHandler(
    filename = function() {
      paste0("sigbenchr_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      shiny::req(results())
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("Install 'rmarkdown' to download the HTML report.", call. = FALSE)
      }
      snap <- results()
      hm_snap <- .snap_method(snap)
      fsig_snap <- .snap_forest_sig(snap)
      m <- report_model()
      man <- reproducibility_manifest()
      summary_unstrat <- if (is.data.frame(m$tier2$summary_tbl)) {
        m$tier2$summary_tbl
      } else {
        suppressWarnings(
          SigBenchR::summarize_auc(snap$res, method = hm_snap, stratify_by = NULL)
        )
      }
      dl_params <- utils::modifyList(
        snap$params,
        list(
          heatmap_method = hm_snap,
          forest_signature = fsig_snap,
          bio_context = m$tier2$ref_axes,
          datasets = snap$datasets,
          signatures = snap$sigs
        )
      )
      SigBenchR::export_report(
        results = snap$res,
        summary_tbl = summary_unstrat,
        sensitivity_tbl = sens(),
        params_list = dl_params,
        output_file = file,
        report_model = m,
        manifest = man
      )
    }
  )

  reproducibility_manifest <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    m <- report_model()
    hm_snap <- .snap_method(snap)
    fsig_snap <- .snap_forest_sig(snap)
    snap_params <- utils::modifyList(
      snap$params,
      list(heatmap_method = hm_snap, forest_signature = fsig_snap)
    )
    SigBenchR::build_reproducibility_manifest(
      results = snap$res,
      summary_tbl = m$tier2$summary_tbl,
      report_model = m,
      params_list = snap_params,
      datasets = snap$datasets,
      signatures = snap$sigs
    )
  })

  output$download_manifest <- shiny::downloadHandler(
    filename = function() {
      paste0("sigbenchr_manifest_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json")
    },
    content = function(file) {
      shiny::req(reproducibility_manifest())
      SigBenchR::export_manifest(reproducibility_manifest(), path = file, format = "json")
    }
  )

  output$download_pub_bundle <- shiny::downloadHandler(
    filename = function() {
      paste0("sigbenchr_publication_bundle_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      shiny::req(results())
      snap <- results()
      hm_snap <- .snap_method(snap)
      fsig_snap <- .snap_forest_sig(snap)
      m <- report_model()
      man <- reproducibility_manifest()
      dl_params <- utils::modifyList(
        snap$params,
        list(
          heatmap_method = hm_snap,
          forest_signature = fsig_snap,
          datasets = snap$datasets,
          signatures = snap$sigs
        )
      )
      SigBenchR::export_publication_bundle(
        results = snap$res,
        output_path = file,
        summary_tbl = m$tier2$summary_tbl,
        report_model = m,
        manifest = man,
        params_list = dl_params,
        method = hm_snap,
        signature = fsig_snap,
        coverage_threshold = snap$params$coverage_threshold,
        zip = TRUE
      )
    }
  )

  .format_summary_display <- function(sm) {
    if (!is.data.frame(sm) || nrow(sm) == 0L) {
      return(sm)
    }
    strata_cols <- intersect(c("treatment", "endpoint_type"), names(sm))
    SigBenchR::format_meta_summary(
      sm,
      stratum_col = if (length(strata_cols) > 0L) strata_cols else NULL
    )
  }

  .translate_analysis_status <- function(st) {
    vapply(as.character(st), function(s) {
      if (is.na(s) || !nzchar(s)) return("Not evaluated (NA)")
      switch(
        s,
        ok = "Valid adjusted model (EPV >= 5, VIF < 5, converged)",
        exploratory_low_epv_or_collinear = "Exploratory: low EPV (< 5) or elevated collinearity (5 <= VIF < 10)",
        failed_collinearity = "Aborted: severe multicollinearity with infiltration axis (VIF >= 10)",
        failed_separation_or_convergence = "Unreliable: quasi-complete separation or GLM non-convergence",
        s
      )
    }, character(1), USE.NAMES = FALSE)
  }

  .translate_background_type <- function(bg) {
    vapply(as.character(bg), function(b) {
      if (is.na(b) || !nzchar(b)) return("Not computed (compute_null = FALSE)")
      switch(
        b,
        subsampled_panel_cohort = "Subsampled panel cohort (~2,800 genes)",
        full_transcriptome = "Full transcriptome (>= 5,000 genes)",
        partial_ortholog_subset = "Partial ortholog-mapped subset (restricted background)",
        b
      )
    }, character(1), USE.NAMES = FALSE)
  }

  .dt_round <- function(df) {
    if (!requireNamespace("DT", quietly = TRUE)) {
      return(df)
    }
    int_like <- c("n", "n_pos", "n_neg", "n_found", "n_total", "n_boot_valid", "universe_size", "n_samples", "n_events", "n_cohorts", "n_datasets")
    num_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], int_like)
    dt <- DT::datatable(
      df,
      options = list(scrollX = TRUE, autoWidth = TRUE, pageLength = 8),
      rownames = FALSE
    )
    if (length(num_cols) > 0L) {
      dt <- DT::formatRound(dt, columns = num_cols, digits = 3)
    }
    dt
  }

  output$overview_stat_alerts <- shiny::renderUI({
    shiny::req(results())
    snap <- results()
    df <- snap$res
    if (!is.data.frame(df) || nrow(df) == 0L) return(NULL)

    items <- list()
    if (isTRUE(snap$params$compute_pr)) {
      if ("small_positive_class" %in% names(df) && any(df$small_positive_class %in% TRUE, na.rm = TRUE)) {
        sc_rows <- unique(df$dataset[df$small_positive_class %in% TRUE])
        items[[length(items) + 1L]] <- shiny::div(
          class = "sb-alert-item",
          shiny::tags$strong("[PR-SMALL-CLASS] Small positive/negative class warning: "),
          sprintf(
            "Cohort(s) %s have fewer than 5 responder or non-responder samples (n_pos < 5 or n_neg < 5). PR-AUC point estimates and stratified bootstrap intervals require cautious interpretation.",
            paste(sc_rows, collapse = ", ")
          )
        )
      }
      if ("ci_unreliable" %in% names(df) && any(df$ci_unreliable %in% TRUE, na.rm = TRUE)) {
        n_unrel <- sum(df$ci_unreliable %in% TRUE, na.rm = TRUE)
        items[[length(items) + 1L]] <- shiny::div(
          class = "sb-alert-item",
          shiny::tags$strong("[PR-CI-UNSTABLE] Unstable PR-AUC confidence interval: "),
          sprintf(
            "%d signature x cohort evaluation(s) flagged ci_unreliable = TRUE due to small minority class size, heavy score ties (>50%%), or <50%% valid bootstrap resamples (n_boot_valid). Do not treat these intervals as confirmatory.",
            n_unrel
          )
        )
      }
    }

    if (isTRUE(snap$params$compute_null) && "background_type" %in% names(df)) {
      bgs <- unique(stats::na.omit(as.character(df$background_type)))
      if (length(bgs) > 0L) {
        bg_labels <- unique(.translate_background_type(bgs))
        is_part <- "partial_ortholog_subset" %in% bgs
        items[[length(items) + 1L]] <- shiny::div(
          class = if (is_part) "sb-alert-item sb-alert-danger" else "sb-alert-item sb-alert-info",
          shiny::tags$strong("[NULL-UNIVERSE] Null permutation background universe: "),
          sprintf(
            "Random gene sets were sampled from: %s. Both within-method BH FDR (null_p_adj_bh) and across-method global BH FDR (null_p_adj_bh_global) are reported below.",
            paste(bg_labels, collapse = "; ")
          )
        )
      }
    }

    if (length(items) == 0L) return(NULL)
    shiny::div(class = "sb-alert-strip", items)
  })

  output$table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      df <- results()$res
      if ("background_type" %in% names(df)) {
        df$background_type_label <- .translate_background_type(df$background_type)
      }
      cols <- c(
        "dataset", "signature", "method", "n", "n_pos", "n_neg",
        "coverage_label", "found_genes", "missing_genes", "low_coverage",
        "auc", "auc_ci_low", "auc_ci_high", "auc_var", "var_zero", "small_cohort",
        "pr_auc", "pr_baseline", "pr_ci_low", "pr_ci_high",
        "small_positive_class", "n_boot_valid", "ci_unreliable",
        "null_p", "null_p_adj_bh", "null_p_adj_bh_global",
        "background_type", "background_type_label",
        "independent", "skipped", "skip_reason"
      )
      .dt_round(df[, intersect(cols, names(df)), drop = FALSE])
    })
  } else NULL

  output$table_fallback <- shiny::renderTable({
    shiny::req(results())
    results()$res[, intersect(
      c(
        "dataset", "signature", "method", "auc", "auc_ci_low", "auc_ci_high",
        "pr_auc", "pr_ci_low", "pr_ci_high", "small_positive_class", "ci_unreliable",
        "null_p", "null_p_adj_bh", "null_p_adj_bh_global", "background_type",
        "low_coverage", "missing_genes"
      ),
      names(results()$res)
    ), drop = FALSE]
  })

  output$coverage_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      df <- results()$res
      cols <- c(
        "dataset", "signature", "method",
        "n_found", "n_total", "coverage_pct", "coverage_label",
        "low_coverage", "found_genes", "missing_genes"
      )
      cov <- df[, intersect(cols, names(df)), drop = FALSE]
      if (all(c("dataset", "signature") %in% names(cov))) {
        cov <- cov[!duplicated(cov[, c("dataset", "signature"), drop = FALSE]), , drop = FALSE]
      }
      .dt_round(cov)
    })
  } else NULL
  output$coverage_fallback <- shiny::renderTable({
    shiny::req(results())
    df <- results()$res
    cols <- intersect(
      c("dataset", "signature", "coverage_label", "low_coverage", "missing_genes"),
      names(df)
    )
    cov <- df[, cols, drop = FALSE]
    if (all(c("dataset", "signature") %in% names(cov))) {
      cov <- cov[!duplicated(cov[, c("dataset", "signature"), drop = FALSE]), , drop = FALSE]
    }
    cov
  })

  output$heatmap <- shiny::renderPlot({
    shiny::req(results())
    snap <- results()
    df <- snap$res
    hm <- .validated_hm_method(snap)
    shiny::validate(shiny::need(
      any(df$method == hm & !(df$skipped %in% TRUE)),
      "No AUC available."
    ))
    SigBenchR::plot_auc_heatmap(
      df,
      method = hm,
      coverage_threshold = snap$params$coverage_threshold
    )
  })

  output$sensitivity_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(sens()); .dt_round(sens()) })
  } else NULL
  output$sensitivity_fallback <- shiny::renderTable({ shiny::req(sens()); sens() })

  output$roc_plot <- shiny::renderPlot({
    if (is.null(input$run) || identical(as.integer(input$run), 0L) || is.null(results())) {
      shiny::validate(shiny::need(FALSE, "Run the analysis before viewing ROC curves."))
    }
    shiny::req(input$roc_dataset)
    snap <- results()
    ds_nm <- as.character(input$roc_dataset)[1]
    hm <- .validated_hm_method(snap)
    shiny::validate(shiny::need(
      !is.na(ds_nm) && nzchar(ds_nm) && ds_nm %in% names(snap$datasets),
      paste0("Selected cohort '", ds_nm, "' is not in the last run snapshot; press 'Run analysis'.")
    ))
    SigBenchR::plot_roc_curves(
      datasets = snap$datasets,
      signatures = snap$sigs,
      dataset_name = ds_nm,
      method = hm,
      direction = snap$params$direction,
      sd_handling = snap$params$sd_handling
    )
  })

  cmp_res <- shiny::eventReactive(input$run_cmp, {
    snap <- results()
    shiny::validate(
      shiny::need(!is.null(snap), "Run the analysis before running DeLong comparison."),
      shiny::need(
        !identical(input$cmp_sig1, input$cmp_sig2),
        "Select two different signatures for DeLong (Signature 1 != Signature 2)."
      )
    )
    shiny::withProgress(message = "DeLong comparison...", value = 0.4, {
      sigs <- snap$sigs
      ds_pool <- snap$datasets
      p <- snap$params

      shiny::validate(
        shiny::need(
          input$cmp_sig1 %in% names(sigs),
          paste0("Signature 1 not in last run snapshot: ", input$cmp_sig1)
        ),
        shiny::need(
          input$cmp_sig2 %in% names(sigs),
          paste0("Signature 2 not in last run snapshot: ", input$cmp_sig2)
        ),
        shiny::need(
          input$cmp_dataset %in% names(ds_pool),
          paste0("Dataset not in last run snapshot: ", input$cmp_dataset)
        )
      )

      sig1 <- sigs[[input$cmp_sig1]]
      sig2 <- sigs[[input$cmp_sig2]]
      g1_set <- sort(unique(toupper(c(sig1$genes, sig1$genes_up, sig1$genes_down))))
      g2_set <- sort(unique(toupper(c(sig2$genes, sig2$genes_up, sig2$genes_down))))
      shiny::validate(
        shiny::need(
          !identical(g1_set, g2_set),
          "Signature 1 and Signature 2 have identical gene sets; self-comparison is not informative."
        )
      )
      ds_use <- ds_pool[[input$cmp_dataset]]
      meth <- if (!is.null(input$heatmap_method) && input$heatmap_method %in% p$methods) {
        input$heatmap_method
      } else {
        p$methods[1]
      }
      dir1 <- SigBenchR::normalize_roc_direction(sig1$direction, default = p$direction)
      dir2 <- SigBenchR::normalize_roc_direction(sig2$direction, default = p$direction)

      s1 <- tryCatch(
        SigBenchR::score_for_dataset(
          ds_use,
          genes = sig1$genes,
          genes_up = sig1$genes_up,
          genes_down = sig1$genes_down,
          method = meth,
          direction = dir1,
          sd_handling = p$sd_handling,
          signature = sig1
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(
        !inherits(s1, "error"),
        paste0("Could not score Signature 1: ", if (inherits(s1, "error")) conditionMessage(s1) else "")
      ))

      s2 <- tryCatch(
        SigBenchR::score_for_dataset(
          ds_use,
          genes = sig2$genes,
          genes_up = sig2$genes_up,
          genes_down = sig2$genes_down,
          method = meth,
          direction = dir2,
          sd_handling = p$sd_handling,
          signature = sig2
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(
        !inherits(s2, "error"),
        paste0("Could not score Signature 2: ", if (inherits(s2, "error")) conditionMessage(s2) else "")
      ))

      cmp <- tryCatch(
        SigBenchR::compare_signatures(
          scores1 = s1$scores,
          scores2 = s2$scores,
          response = s1$response,
          dataset1 = input$cmp_dataset,
          dataset2 = input$cmp_dataset,
          direction1 = dir1,
          direction2 = dir2,
          name1 = input$cmp_sig1,
          name2 = input$cmp_sig2
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(
        !inherits(cmp, "error"),
        paste0("DeLong comparison error: ", if (inherits(cmp, "error")) conditionMessage(cmp) else "")
      ))
      cmp
    })
  }, ignoreInit = TRUE)

  output$cmp_status <- shiny::renderText({
    if (identical(input$cmp_sig1, input$cmp_sig2)) {
      return("Warning: Signature 1 and Signature 2 must differ.")
    }
    if (is.null(input$run_cmp) || identical(as.integer(input$run_cmp), 0L) || is.null(cmp_res())) {
      return("Select two signatures on the same dataset and run DeLong.")
    }
    r <- cmp_res()
    sprintf(
      "p=%.4g | p_adj_bh=%.4g (equals raw p for a single comparison) | AUC1=%.3f | AUC2=%.3f",
      r$p_value, r$p_adj_bh, r$auc1, r$auc2
    )
  })
  output$cmp_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(cmp_res()); .dt_round(cmp_res()) })
  } else NULL
  output$cmp_fallback <- shiny::renderTable({ shiny::req(cmp_res()); cmp_res() })

  output$summary_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(summary_res())
      .dt_round(.format_summary_display(summary_res()))
    })
  } else NULL
  output$summary_fallback <- shiny::renderTable({
    shiny::req(summary_res())
    .format_summary_display(summary_res())
  })

  output$forest_plot <- shiny::renderPlot({
    shiny::req(results())
    snap <- results()
    hm <- .validated_hm_method(snap)
    fsig <- if (!is.null(input$forest_sig) &&
                length(input$forest_sig) >= 1L &&
                nzchar(as.character(input$forest_sig)[1])) {
      as.character(input$forest_sig)[1]
    } else {
      .snap_forest_sig(snap)
    }
    shiny::validate(
      shiny::need(
        !is.na(fsig) && nzchar(fsig) && fsig %in% names(snap$sigs),
        paste0("Selected forest signature '", fsig, "' is not in the last run snapshot; press 'Run analysis'.")
      ),
      shiny::need(
        any(snap$res$signature == fsig &
              snap$res$method == hm &
              !(snap$res$skipped %in% TRUE)),
        "No AUC for the selected signature/method."
      )
    )
    SigBenchR::plot_auc_forest(
      snap$res,
      signature = fsig,
      method = hm,
      summary = summary_res_forest()
    )
  })

  output$jaccard_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      snap <- results()
      shiny::validate(shiny::need(
        length(snap$sigs) >= 2L,
        "Select at least two signatures for Jaccard overlap."
      ))
      .dt_round(SigBenchR::jaccard_signatures(snap$sigs))
    })
  } else NULL
  output$jaccard_fallback <- shiny::renderTable({
    shiny::req(results())
    snap <- results()
    shiny::validate(shiny::need(
      length(snap$sigs) >= 2L,
      "Select at least two signatures for Jaccard overlap."
    ))
    SigBenchR::jaccard_signatures(snap$sigs)
  })

  output$corr_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      snap <- results()
      hm <- .validated_hm_method(snap)
      shiny::validate(shiny::need(
        length(snap$sigs) >= 2L &&
          length(unique(snap$res$signature[!(snap$res$skipped %in% TRUE)])) >= 2L,
        "Run at least two signatures for correlation."
      ))
      .dt_round(SigBenchR::auc_correlation(snap$res, method = hm))
    })
  } else NULL
  output$corr_fallback <- shiny::renderTable({
    shiny::req(results())
    snap <- results()
    hm <- .validated_hm_method(snap)
    shiny::validate(shiny::need(
      length(snap$sigs) >= 2L &&
        length(unique(snap$res$signature[!(snap$res$skipped %in% TRUE)])) >= 2L,
      "Run at least two signatures for correlation."
    ))
    SigBenchR::auc_correlation(snap$res, method = hm)
  })

  gene_contrib_res <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    sig_nm <- if (!is.null(input$gene_contrib_sig) &&
                  length(input$gene_contrib_sig) >= 1L &&
                  nzchar(as.character(input$gene_contrib_sig)[1])) {
      as.character(input$gene_contrib_sig)[1]
    } else {
      names(snap$sigs)[1]
    }
    shiny::validate(shiny::need(
      sig_nm %in% names(snap$sigs),
      paste0("Selected signature '", sig_nm, "' is not in the last run snapshot; press 'Run analysis'.")
    ))
    hm_active <- .validated_hm_method(snap)
    SigBenchR::analyze_gene_contributions(
      datasets = snap$datasets,
      signature = snap$sigs[[sig_nm]],
      method = hm_active,
      sd_handling = snap$params$sd_handling,
      direction = snap$params$direction
    )
  })

  output$gene_contrib_plot <- shiny::renderPlot({
    gc <- gene_contrib_res()
    shiny::validate(shiny::need(
      is.data.frame(gc$summary) && nrow(gc$summary) > 0L,
      "Not enough matched genes (>= 2) in the selected signature/snapshot."
    ))
    SigBenchR::plot_gene_contributions(gc)
  })

  output$gene_contrib_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      gc <- gene_contrib_res()
      .dt_round(gc$summary)
    })
  } else NULL
  output$gene_contrib_fallback <- shiny::renderTable({
    gene_contrib_res()$summary
  })

  surv_bench_res <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    hm <- .validated_hm_method(snap)
    suppressWarnings(
      tryCatch(
        SigBenchR::run_survival_benchmark(
          datasets = snap$datasets,
          signatures = snap$sigs,
          method = hm
        ),
        error = function(e) data.frame()
      )
    )
  })

  output$survival_alerts <- shiny::renderUI({
    shiny::req(results())
    sb <- surv_bench_res()
    if (!is.data.frame(sb) || nrow(sb) == 0L) {
      return(shiny::div(
        class = "sb-alert-strip",
        shiny::div(
          class = "sb-alert-item sb-alert-info",
          shiny::tags$strong("[SURVIVAL-DATA] No time-to-event columns in active cohorts: "),
          "None of the currently selected snapshot cohorts contain overall survival columns (os_time / os_event). Load a cohort with survival follow-up to populate Cox HR (p_adj_bh) and Log-rank (logrank_p_adj_bh) evaluations."
        )
      ))
    }
    items <- list()
    if ("ph_violated" %in% names(sb) && any(sb$ph_violated %in% TRUE, na.rm = TRUE)) {
      n_ph <- sum(sb$ph_violated %in% TRUE, na.rm = TRUE)
      items[[length(items) + 1L]] <- shiny::div(
        class = "sb-alert-item",
        shiny::tags$strong("[SURV-PH-DIAGNOSTIC] Proportional hazards assumption warning: "),
        sprintf(
          "%d evaluation(s) have unadjusted Schoenfeld residual test ph_p < 0.05 (ph_violated = TRUE). Interpret constant Cox hazard ratios cautiously; ph_p is left unadjusted intentionally so model violations are not masked by FDR correction.",
          n_ph
        )
      )
    }
    if (length(items) == 0L) return(NULL)
    shiny::div(class = "sb-alert-strip", items)
  })

  output$survival_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      sb <- surv_bench_res()
      cols <- c(
        "cohort", "signature", "method", "endpoint", "n_samples", "n_events",
        "hr", "hr_ci_low", "hr_ci_high", "cox_p", "p_adj_bh", "c_index",
        "logrank_p", "logrank_p_adj_bh", "median_surv_high", "median_surv_low",
        "ph_p", "ph_violated", "small_sample", "surv_note"
      )
      .dt_round(sb[, intersect(cols, names(sb)), drop = FALSE])
    })
  } else NULL

  output$survival_fallback <- shiny::renderTable({
    shiny::req(results())
    sb <- surv_bench_res()
    cols <- c(
      "cohort", "signature", "method", "n_samples", "n_events",
      "hr", "cox_p", "p_adj_bh", "logrank_p", "logrank_p_adj_bh",
      "ph_p", "ph_violated"
    )
    sb[, intersect(cols, names(sb)), drop = FALSE]
  })

  beyond_inf_res <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    sig_nm <- if (!is.null(input$beyond_inf_sig) &&
                  length(input$beyond_inf_sig) >= 1L &&
                  nzchar(as.character(input$beyond_inf_sig)[1])) {
      as.character(input$beyond_inf_sig)[1]
    } else {
      names(snap$sigs)[1]
    }
    shiny::validate(shiny::need(
      sig_nm %in% names(snap$sigs),
      paste0("Selected signature '", sig_nm, "' is not in the last run snapshot; press 'Run analysis'.")
    ))
    ax <- if (!is.null(input$beyond_inf_axis) && nzchar(input$beyond_inf_axis)) {
      input$beyond_inf_axis
    } else {
      "CD8_T_Cell"
    }
    hm_active <- .validated_hm_method(snap)
    SigBenchR::test_beyond_infiltration(
      datasets = snap$datasets,
      signature = snap$sigs[[sig_nm]],
      infiltration_axis = ax,
      method = hm_active,
      sd_handling = snap$params$sd_handling,
      direction = snap$params$direction
    )
  })

  output$beyond_inf_alerts <- shiny::renderUI({
    bi_all <- beyond_inf_res()
    pc <- bi_all$per_cohort
    sm <- bi_all$summary
    if (!is.data.frame(pc) || nrow(pc) == 0L) return(NULL)

    items <- list()
    if (is.data.frame(sm) && nrow(sm) > 0L && isTRUE(sm$direction_conflict[1])) {
      items[[length(items) + 1L]] <- shiny::div(
        class = "sb-alert-item sb-alert-danger",
        shiny::tags$strong("[DIRECTION-CONFLICT] Cross-cohort effect direction conflict (direction_conflict = TRUE): "),
        sprintf(
          "The adjusted signature effect (beta_adjusted) has opposite signs across cohorts. Directional Stouffer Z combined_lrt_p = %.4g properly accounts for cancellation, whereas unsigned Fisher omnibus p (%.4g) would misleadingly pool opposing effects.",
          sm$combined_lrt_p[1], sm$fisher_omnibus_p[1]
        )
      )
    }

    if ("low_epv" %in% names(pc) && any(pc$low_epv %in% TRUE, na.rm = TRUE)) {
      low_c <- pc$cohort[pc$low_epv %in% TRUE]
      min_epv <- min(pc$epv[pc$low_epv %in% TRUE], na.rm = TRUE)
      items[[length(items) + 1L]] <- shiny::div(
        class = "sb-alert-item",
        shiny::tags$strong("[LOW-EPV] Small sample events-per-variable warning (EPV < 5): "),
        sprintf(
          "Cohort(s) %s have low Events Per Variable (minimum EPV = %.1f < 5). Adjusted logistic regression coefficients and LRT p-values have limited statistical reliability and should be treated as exploratory.",
          paste(low_c, collapse = ", "), min_epv
        )
      )
    }

    if ("high_collinearity" %in% names(pc) && any(pc$high_collinearity %in% TRUE, na.rm = TRUE)) {
      col_c <- pc$cohort[pc$high_collinearity %in% TRUE]
      max_vif <- max(pc$vif[pc$high_collinearity %in% TRUE], na.rm = TRUE)
      aborted <- any(pc$analysis_status == "failed_collinearity", na.rm = TRUE)
      items[[length(items) + 1L]] <- shiny::div(
        class = if (aborted) "sb-alert-item sb-alert-danger" else "sb-alert-item",
        shiny::tags$strong("[HIGH-VIF] Multicollinearity with reference infiltration axis: "),
        sprintf(
          "Cohort(s) %s exhibit high Variance Inflation Factor (max VIF = %.2f >= 5) with %s. %s",
          paste(col_c, collapse = ", "),
          max_vif,
          sm$infiltration_axis[1],
          if (aborted) "Cohorts with VIF >= 10 were aborted (failed_collinearity) because the signature is statistically indistinguishable from baseline infiltration." else "Elevated collinearity inflates beta_adjusted standard errors."
        )
      )
    }

    if (length(items) == 0L) return(NULL)
    shiny::div(class = "sb-alert-strip", items)
  })

  output$beyond_inf_summary <- shiny::renderText({
    bi <- beyond_inf_res()$summary
    if (!is.data.frame(bi) || nrow(bi) == 0L || bi$n_cohorts[1] == 0L) {
      return("No cohorts with sufficient infiltration & signature genes in snapshot.")
    }
    sprintf(
      "Axis: %s | Cohorts: %d | Independent value (BH LRT p_adj < 0.05 & delta > 0): %d/%d | Directional Stouffer p: %.4g (Z = %.2f) | Direction conflict: %s | Mean CV delta-AUC: %.3f",
      bi$infiltration_axis[1], bi$n_cohorts[1], bi$n_independent[1], bi$n_cohorts[1],
      bi$combined_lrt_p[1],
      if ("stouffer_z" %in% names(bi)) bi$stouffer_z[1] else NA_real_,
      if ("direction_conflict" %in% names(bi) && isTRUE(bi$direction_conflict[1])) "YES (opposite signs across cohorts)" else "No",
      bi$mean_delta_auc[1]
    )
  })

  output$beyond_inf_plot <- shiny::renderPlot({
    bi <- beyond_inf_res()$per_cohort
    shiny::validate(shiny::need(
      is.data.frame(bi) && nrow(bi) > 0L && any(is.finite(bi$delta_auc_combined)),
      "No finite cross-validated delta-AUC values available to plot."
    ))
    df <- bi[is.finite(bi$delta_auc_combined), , drop = FALSE]
    df$cohort_label <- sprintf("%s (n=%d, EPV=%.1f)", df$cohort, df$n, df$epv)
    df$cohort_label <- factor(df$cohort_label, levels = rev(unique(df$cohort_label)))
    df$status_label <- .translate_analysis_status(df$analysis_status)

    ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$delta_auc_combined, y = .data$cohort_label, color = .data$status_label)
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$delta_auc_ci_low, xmax = .data$delta_auc_ci_high),
        orientation = "y",
        width = 0.2,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(size = 2.8) +
      ggplot2::labs(
        title = "Paired Cross-Validated Delta-AUC Beyond Infiltration",
        subtitle = "Delta-AUC = AUC(signature + infiltration) - AUC(infiltration alone); error bars: fold-level 95% CI",
        x = "Paired CV Delta-AUC (95% CI)",
        y = "Cohort (sample size & EPV)",
        color = "Model diagnostic status"
      ) +
      ggplot2::theme_bw(base_size = 12)
  })

  .format_beyond_inf_df <- function(pc) {
    if (!is.data.frame(pc) || nrow(pc) == 0L) return(pc)
    pc$analysis_status_explanation <- .translate_analysis_status(pc$analysis_status)
    cols <- c(
      "cohort", "n", "n_pos", "n_neg", "epv", "low_epv",
      "r_with_infiltration", "vif", "high_collinearity", "jaccard_overlap",
      "auc_infiltration_only", "auc_combined",
      "delta_auc_combined", "sd_delta_auc", "delta_auc_ci_low", "delta_auc_ci_high",
      "beta_adjusted", "or_adjusted", "p_adjusted", "lrt_p", "lrt_p_adj_bh",
      "adds_independent_value", "analysis_status", "analysis_status_explanation",
      "converged", "separation_flag"
    )
    pc[, intersect(cols, names(pc)), drop = FALSE]
  }

  output$beyond_inf_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      bi <- beyond_inf_res()
      .dt_round(.format_beyond_inf_df(bi$per_cohort))
    })
  } else NULL
  output$beyond_inf_fallback <- shiny::renderTable({
    .format_beyond_inf_df(beyond_inf_res()$per_cohort)
  })
}

shiny::shinyApp(ui, server)
