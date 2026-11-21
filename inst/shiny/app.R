# SigBenchR Shiny - tek dosya

# Global: oturum basina yeniden yuklemeyi onlemek icin server disinda
data(datasets_ici, package = "SigBenchR", envir = environment())
data(signatures_demo, package = "SigBenchR", envir = environment())

# Fail-closed: yalniz run_app(allow_upload = TRUE) ile acilir
.sigbenchr_allow_upload <- isTRUE(
  shiny::getShinyOption("sigbenchr_allow_upload", default = FALSE)
)

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
      shiny::span(class = "sb-badge", "v0.4.3"),
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
        )
      ),
      shiny::tags$details(
        class = "sb-acc",
        open = NA,
        shiny::tags$summary("1. Signatures and cohorts"),
        shiny::checkboxGroupInput(
          "sig_choice", "Built-in signatures",
          choices = c("Cytokine-5", "IFNG6", "Random3"),
          selected = c("Cytokine-5", "IFNG6", "Random3")
        ),
        shiny::textAreaInput(
          "custom_genes",
          "Custom signature (optional)",
          rows = 3,
          placeholder = "Separate genes by comma, space, or newline"
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
        }
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
            shiny::tags$h4(class = "sb-section", "Results table"),
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
      output$upload_status <- shiny::renderText({
        paste0("Error: ", conditionMessage(ds))
      })
      return()
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
    output$upload_status <- shiny::renderText({
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
      msg
    })
  })

  active_sig_names <- shiny::reactive({
    chosen <- as.character(input$sig_choice)
    chosen <- chosen[nzchar(chosen)]
    custom_txt <- if (is.null(input$custom_genes)) "" else input$custom_genes
    if (nzchar(trimws(custom_txt))) {
      chosen <- c(chosen, "Custom")
    }
    unique(chosen)
  })

  build_sigs <- function(direction, sig_choice, custom_genes) {
    chosen <- as.character(sig_choice)
    chosen <- chosen[nzchar(chosen)]
    sigs <- list()
    for (nm in chosen) {
      shiny::validate(shiny::need(nm %in% names(signatures_demo), paste0("Unknown signature: ", nm)))
      sig <- signatures_demo[[nm]]
      sig$direction <- direction
      sigs[[nm]] <- sig
    }
    custom_txt <- if (is.null(custom_genes)) "" else custom_genes
    if (nzchar(trimws(custom_txt))) {
      genes <- SigBenchR::parse_genes(custom_txt)
      shiny::validate(shiny::need(length(genes) >= 1, "Enter at least one gene for the custom signature."))
      # derived_from omitted -> independent = NA (fail-closed)
      sigs[["Custom"]] <- list(genes = genes, direction = direction)
    }
    shiny::validate(shiny::need(length(sigs) >= 1, "Select at least one built-in signature or add a custom signature."))
    sigs
  }

  current_params <- shiny::reactive({
    list(
      methods = as.character(input$methods),
      direction = as.character(input$direction),
      coverage_threshold = as.numeric(input$coverage_threshold),
      sensitivity_threshold = as.numeric(input$sensitivity_threshold),
      n_boot = as.integer(input$n_boot),
      n_perm = as.integer(input$n_perm),
      seed = as.integer(input$seed),
      compute_pr = isTRUE(input$compute_pr),
      compute_null = isTRUE(input$compute_null),
      sd_handling = as.character(input$sd_handling),
      match_expression = isTRUE(input$match_expression),
      datasets = as.character(input$datasets),
      sig_choice = as.character(input$sig_choice),
      custom_genes = if (is.null(input$custom_genes)) "" else as.character(input$custom_genes),
      heatmap_method = as.character(input$heatmap_method),
      forest_signature = as.character(input$forest_sig)
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
      sigs <- c("Cytokine-5", "IFNG6", "Random3")
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
    sigs <- build_sigs(params$direction, params$sig_choice, params$custom_genes)
    pool <- active_dataset_pool()
    miss_ds <- setdiff(params$datasets, names(pool))
    shiny::validate(shiny::need(
      length(miss_ds) == 0L,
      paste0("Unknown dataset: ", paste(miss_ds, collapse = ", "))
    ))
    ds <- pool[params$datasets]

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
      "methods", "direction", "coverage_threshold", "n_boot", "n_perm", "seed",
      "compute_pr", "compute_null", "sd_handling", "match_expression",
      "datasets", "sig_choice", "custom_genes"
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

    df <- snap$res
    p <- snap$params
    n_ds <- length(unique(as.character(p$datasets)))
    n_sig <- length(snap$sigs)
    hm <- if (!is.null(input$heatmap_method) && nzchar(input$heatmap_method)) {
      input$heatmap_method
    } else {
      p$heatmap_method
    }
    sd_lab <- switch(
      as.character(p$sd_handling),
      exclude = "exclude",
      as_0 = "as_0 (ORR)",
      as_1 = "as_1 (DCR)",
      as.character(p$sd_handling)
    )

    indep_flag <- as.logical(df$independent)
    indep <- df[
      df$method == hm &
        !(df$skipped %in% TRUE) &
        !is.na(indep_flag) & indep_flag &
        !(df$low_coverage %in% TRUE) &
        is.finite(df$auc),
      ,
      drop = FALSE
    ]
    best_auc <- if (nrow(indep) > 0L) {
      i <- which.max(indep$auc)
      sprintf("%.3f (%s · %s)", indep$auc[i], indep$signature[i], indep$dataset[i])
    } else {
      "No independent result"
    }

    n_low <- sum(df$low_coverage %in% TRUE, na.rm = TRUE)
    thr_pct <- round(100 * as.numeric(p$coverage_threshold))

    shiny::div(
      class = "sb-kpi-row",
      shiny::div(
        class = "sb-kpi",
        shiny::tags$p(class = "sb-kpi-label", "Cohort / Signature"),
        shiny::tags$p(class = "sb-kpi-value", sprintf("%d cohorts · %d signatures", n_ds, n_sig))
      ),
      shiny::div(
        class = "sb-kpi",
        shiny::tags$p(class = "sb-kpi-label", "Method & SD"),
        shiny::tags$p(class = "sb-kpi-value", sprintf("%s · %s", hm, sd_lab))
      ),
      shiny::div(
        class = "sb-kpi",
        shiny::tags$p(class = "sb-kpi-label", "Best independent AUC"),
        shiny::tags$p(class = "sb-kpi-value", best_auc)
      ),
      shiny::div(
        class = "sb-kpi",
        shiny::tags$p(class = "sb-kpi-label", sprintf("Low coverage (< %d%%)", thr_pct)),
        shiny::tags$p(class = "sb-kpi-value", sprintf("%d cells", as.integer(n_low)))
      )
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

  summary_res <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    hm <- if (!is.null(input$heatmap_method) && nzchar(input$heatmap_method)) {
      input$heatmap_method
    } else {
      snap$params$heatmap_method
    }
    strat <- if (isTRUE(input$stratify_by)) {
      c("treatment", "endpoint_type")
    } else {
      NULL
    }
    suppressWarnings(
      SigBenchR::summarize_auc(snap$res, method = hm, stratify_by = strat)
    )
  })

  # Unstratified summary for forest reference (avoids multi-row stratify_by clash)
  summary_res_forest <- shiny::reactive({
    shiny::req(results())
    snap <- results()
    hm <- if (!is.null(input$heatmap_method) && nzchar(input$heatmap_method)) {
      input$heatmap_method
    } else {
      snap$params$heatmap_method
    }
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
    custom_txt <- if (is.null(p$custom_genes)) "" else p$custom_genes
    genes <- if (nzchar(trimws(custom_txt))) {
      SigBenchR::parse_genes(custom_txt)
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
        match_expression = p$match_expression
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
      hm_active <- if (!is.null(input$heatmap_method) && nzchar(input$heatmap_method)) {
        input$heatmap_method
      } else {
        snap$params$heatmap_method
      }
      summary_unstrat <- suppressWarnings(
        SigBenchR::summarize_auc(snap$res, method = hm_active, stratify_by = NULL)
      )
      bio_ctx <- tryCatch(
        SigBenchR::correlate_reference_axes(
          datasets = snap$datasets,
          signatures = snap$sigs,
          method = hm_active,
          sd_handling = snap$params$sd_handling
        ),
        error = function(e) NULL
      )
      SigBenchR::export_report(
        results = snap$res,
        summary_tbl = summary_unstrat,
        sensitivity_tbl = sens(),
        params_list = list(
          methods = snap$params$methods,
          direction = snap$params$direction,
          coverage_threshold = snap$params$coverage_threshold,
          sensitivity_threshold = snap$params$sensitivity_threshold,
          n_boot = snap$params$n_boot,
          n_perm = snap$params$n_perm,
          seed = snap$params$seed,
          compute_null = isTRUE(snap$params$compute_null),
          sd_handling = snap$params$sd_handling,
          match_expression = snap$params$match_expression,
          heatmap_method = hm_active,
          forest_signature = snap$params$forest_signature,
          bio_context = bio_ctx
        ),
        output_file = file
      )
    }
  )

  .em <- "\u2014"
  .en <- "\u2013"

  .fmt3 <- function(x) {
    x <- as.numeric(x)
    ifelse(is.finite(x), sprintf("%.3f", x), .em)
  }

  .fmt_i2 <- function(x) {
    x <- as.numeric(x)
    ifelse(is.finite(x), sprintf("%.1f", x), .em)
  }

  .fmt_auc_ci <- function(est, lo, hi) {
    est <- as.numeric(est)
    lo <- as.numeric(lo)
    hi <- as.numeric(hi)
    out <- rep(.em, length(est))
    full <- is.finite(est) & is.finite(lo) & is.finite(hi)
    out[full] <- sprintf("%.3f (%.3f%s%.3f)", est[full], lo[full], .en, hi[full])
    point_only <- is.finite(est) & !full
    out[point_only] <- sprintf("%.3f", est[point_only])
    out
  }

  .translate_note <- function(x) {
    x <- as.character(x)
    miss <- is.na(x) | !nzchar(x)
    x[miss] <- .em
    x <- gsub(
      "([0-9]+) bagimsiz-olmayan/NA satir ozetten cikarildi",
      "\\1 non-independent/NA cohorts excluded",
      x
    )
    x <- gsub(
      "([0-9]+) dusuk kapsama satir ozetten cikarildi",
      "\\1 low-coverage cohorts excluded",
      x
    )
    x <- gsub(
      "([0-9]+) kucuk kohort uyarisi",
      "\\1 small-cohort warning(s)",
      x
    )
    x <- gsub(
      "([0-9]+) var_zero/NA varyans satir IVW'den cikarildi",
      "\\1 zero/NA-variance rows excluded from IVW",
      x
    )
    x
  }

  .format_summary_display <- function(sm) {
    if (!is.data.frame(sm) || nrow(sm) == 0L) {
      return(sm)
    }
    k_col <- if ("n_cohorts" %in% names(sm)) sm$n_cohorts else sm$n_datasets
    has_strata <- any(c("treatment", "endpoint_type") %in% names(sm))
    if (isTRUE(has_strata)) {
      stratum <- if ("treatment" %in% names(sm) && "endpoint_type" %in% names(sm)) {
        paste(as.character(sm$treatment), as.character(sm$endpoint_type), sep = " / ")
      } else if ("treatment" %in% names(sm)) {
        as.character(sm$treatment)
      } else {
        as.character(sm$endpoint_type)
      }
      out <- data.frame(
        Signature = as.character(sm$signature),
        Method = as.character(sm$method),
        Stratum = stratum,
        k = as.integer(k_col),
        `N-Weighted AUC` = .fmt3(sm$auc_n_weighted),
        `REML AUC (95% CI)` = .fmt_auc_ci(sm$auc_re, sm$auc_re_ci_low, sm$auc_re_ci_high),
        `I2 (%)` = .fmt_i2(sm$I2),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      out <- data.frame(
        Signature = as.character(sm$signature),
        Method = as.character(sm$method),
        k = as.integer(k_col),
        `N-Weighted AUC` = .fmt3(sm$auc_n_weighted),
        `IVW AUC` = .fmt3(sm$auc_ivw),
        `REML AUC (95% CI)` = .fmt_auc_ci(sm$auc_re, sm$auc_re_ci_low, sm$auc_re_ci_high),
        `I2 (%)` = .fmt_i2(sm$I2),
        Note = if ("note" %in% names(sm)) .translate_note(sm$note) else rep(.em, nrow(sm)),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
    out
  }

  .dt_round <- function(df) {
    if (!requireNamespace("DT", quietly = TRUE)) {
      return(df)
    }
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
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

  output$table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      df <- results()$res
      cols <- c(
        "dataset", "signature", "method", "n", "n_pos", "n_neg",
        "coverage_label", "found_genes", "missing_genes", "low_coverage",
        "auc", "auc_ci_low", "auc_ci_high", "auc_var", "var_zero", "small_cohort",
        "pr_auc", "pr_baseline", "pr_ci_low", "pr_ci_high",
        "null_p", "null_p_adj_bh", "independent",
        "skipped", "skip_reason"
      )
      .dt_round(df[, intersect(cols, names(df)), drop = FALSE])
    })
  } else NULL

  output$table_fallback <- shiny::renderTable({
    shiny::req(results())
    results()$res[, intersect(
      c("dataset", "signature", "method", "auc", "pr_auc", "null_p", "low_coverage", "missing_genes"),
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
    hm <- input$heatmap_method
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
    SigBenchR::plot_roc_curves(
      datasets = snap$datasets,
      signatures = snap$sigs,
      dataset_name = input$roc_dataset,
      method = input$heatmap_method,
      direction = snap$params$direction,
      sd_handling = snap$params$sd_handling
    )
  })

  cmp_res <- shiny::eventReactive(input$run_cmp, {
    shiny::validate(
      shiny::need(
        !identical(input$cmp_sig1, input$cmp_sig2),
        "Select two different signatures for DeLong (Signature 1 != Signature 2)."
      )
    )
    shiny::withProgress(message = "DeLong comparison...", value = 0.4, {
      genes1 <- if (identical(input$cmp_sig1, "Custom")) {
        SigBenchR::parse_genes(input$custom_genes)
      } else {
        signatures_demo[[input$cmp_sig1]]$genes
      }
      genes2 <- if (identical(input$cmp_sig2, "Custom")) {
        SigBenchR::parse_genes(input$custom_genes)
      } else {
        signatures_demo[[input$cmp_sig2]]$genes
      }
      shiny::validate(
        shiny::need(length(genes1) >= 1, "Signature 1 genes not found."),
        shiny::need(length(genes2) >= 1, "Signature 2 genes not found.")
      )

      cmp_pool <- active_dataset_pool()
      shiny::validate(shiny::need(
        input$cmp_dataset %in% names(cmp_pool),
        paste0("Dataset not found: ", input$cmp_dataset)
      ))
      s1 <- tryCatch(
        SigBenchR::score_for_dataset(
          cmp_pool[[input$cmp_dataset]],
          genes = genes1,
          method = input$heatmap_method,
          direction = input$direction,
          sd_handling = input$sd_handling
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(
        !inherits(s1, "error"),
        paste0("Could not score Signature 1: ", if (inherits(s1, "error")) conditionMessage(s1) else "")
      ))

      s2 <- tryCatch(
        SigBenchR::score_for_dataset(
          cmp_pool[[input$cmp_dataset]],
          genes = genes2,
          method = input$heatmap_method,
          direction = input$direction,
          sd_handling = input$sd_handling
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
          direction1 = input$direction,
          direction2 = input$direction,
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
    shiny::validate(shiny::need(
      any(snap$res$signature == input$forest_sig &
            snap$res$method == input$heatmap_method &
            !(snap$res$skipped %in% TRUE)),
      "No AUC for the selected signature/method."
    ))
    SigBenchR::plot_auc_forest(
      snap$res,
      signature = input$forest_sig,
      method = input$heatmap_method,
      summary = summary_res_forest()
    )
  })

  output$jaccard_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      .dt_round(SigBenchR::jaccard_signatures(results()$sigs))
    })
  } else NULL
  output$jaccard_fallback <- shiny::renderTable({
    shiny::req(results())
    SigBenchR::jaccard_signatures(results()$sigs)
  })

  output$corr_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      snap <- results()
      shiny::validate(shiny::need(
        length(unique(snap$res$signature[!(snap$res$skipped %in% TRUE)])) >= 2,
        "Run at least two signatures for correlation."
      ))
      .dt_round(SigBenchR::auc_correlation(snap$res, method = input$heatmap_method))
    })
  } else NULL
  output$corr_fallback <- shiny::renderTable({
    shiny::req(results())
    snap <- results()
    SigBenchR::auc_correlation(snap$res, method = input$heatmap_method)
  })
}

shiny::shinyApp(ui, server)
