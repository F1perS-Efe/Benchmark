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
    ))
  ),
  shiny::titlePanel("SigBenchR — Gen İmzası Doğrulama"),
  shiny::tags$p(
    shiny::em(
      "Dahili halka açık melanoma ICI kohortları üzerinde keşifsel kıyaslama aracı — yalnızca gen sembolleri girin. ",
      "Araştırma amaçlıdır; klinik karar aracı değildir."
    )
  ),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::checkboxGroupInput(
        "sig_choice", "Hazır imzalar",
        choices = c("Cytokine-5", "IFNG6", "Random3"),
        selected = c("Cytokine-5", "IFNG6", "Random3")
      ),
      shiny::textAreaInput(
        "custom_genes",
        "Özel imza ekle — opsiyonel",
        rows = 3,
        placeholder = "Genleri virgül, boşluk veya satırla ayırın"
      ),
      shiny::checkboxGroupInput(
        "methods", "Skorlama yöntemleri",
        choices = c("mean", "zscore", "median", "rank"),
        selected = c("mean", "zscore")
      ),
      shiny::radioButtons(
        "heatmap_method", "Aktif Görünüm / Özet Yöntemi",
        choices = c("mean", "zscore"), selected = "mean"
      ),
      shiny::radioButtons(
        "direction", "İmza yönü",
        choices = c("Yüksek=yanıtlı (<)" = "<", "Düşük=yanıtlı (>)" = ">"),
        selected = "<"
      ),
      shiny::radioButtons(
        "sd_handling", "SD politikası",
        choices = c(
          "SD dışla [CR/PR vs PD]" = "exclude",
          "ORR [CR/PR vs SD/PD]" = "as_0",
          "DCR [CR/PR/SD vs PD]" = "as_1"
        ),
        selected = "exclude"
      ),
      shiny::checkboxInput(
        "match_expression",
        "Null: ifade düzeyine eşleştirilmiş örnekleme",
        FALSE
      ),
      shiny::sliderInput("coverage_threshold", "Kapsama eşiği", 0.5, 1, 0.70, 0.05),
      shiny::sliderInput("sensitivity_threshold", "Dayanıklılık eşiği (keyfi)", 0.05, 0.3, 0.10, 0.01),
      shiny::checkboxInput("compute_pr", "PR-AUC + bootstrap GA", TRUE),
      shiny::numericInput("n_boot", "PR bootstrap (n_boot)", value = 200, min = 0, max = 2000, step = 50),
      shiny::checkboxInput("compute_null", "Null permütasyon", FALSE),
      shiny::numericInput("n_perm", "Null B", value = 1000, min = 10, max = 2000, step = 50),
      shiny::numericInput("seed", "Seed", value = 42, min = 1, step = 1),
      shiny::checkboxGroupInput(
        "datasets", "Veri setleri (ICI kohortları)",
        choices = c(
          "van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019"
        ),
        selected = c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
      ),
      if (isTRUE(.sigbenchr_allow_upload)) {
        shiny::tagList(
          shiny::hr(),
          shiny::h5("Yerel Kohort Ekle (CSV/TSV)"),
          shiny::fileInput("upload_expr", "Ifade matrisi (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
          shiny::fileInput("upload_meta", "Klinik metadata (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
          shiny::checkboxInput("upload_genes_in_rows", "Genler satırlarda", TRUE),
          shiny::selectInput(
            "upload_gene_id_type", "Gen ID tipi",
            choices = c("HGNC", "Ensembl", "other"), selected = "HGNC"
          ),
          shiny::selectInput(
            "upload_normalization", "Normalizasyon",
            choices = c("TPM", "FPKM", "normalized_counts", "log_normalized", "other"),
            selected = "TPM"
          ),
          shiny::selectInput(
            "upload_scale", "Ölçek",
            choices = c("linear", "log2"), selected = "log2"
          ),
          shiny::textInput("upload_sample_id_col", "sample_id sütunu", value = "sample_id"),
          shiny::textInput("upload_response_col", "response sütunu", value = "response_label"),
          shiny::textInput("upload_dataset_name", "Kohort adı", value = "user_cohort"),
          shiny::actionButton("upload_add", "Kohortu Doğrula ve Ekle", class = "btn-info"),
          shiny::verbatimTextOutput("upload_status")
        )
      } else {
        shiny::tags$p(
          shiny::em(
            "Özel kohort yükleme hasta gizliliği için yalnızca yerel R oturumunda (`SigBenchR::run_app(allow_upload = TRUE)`) aktiftir."
          )
        )
      },
      shiny::actionButton("run", "Analizi çalıştır", class = "btn-primary"),
      shiny::downloadButton("download_report", "Tüm sonuçları HTML rapor olarak indir"),
      shiny::hr(),
      shiny::h5("DeLong: iki imza (aynı veri seti)"),
      shiny::selectInput(
        "cmp_dataset", "Veri seti",
        choices = c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
      ),
      shiny::selectInput("cmp_sig1", "İmza 1", choices = c("Cytokine-5", "IFNG6", "Random3")),
      shiny::selectInput("cmp_sig2", "İmza 2", choices = c("IFNG6", "Cytokine-5", "Random3")),
      shiny::actionButton("run_cmp", "DeLong karşılaştır")
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Sonuçlar",
          shiny::verbatimTextOutput("status"),
          shiny::uiOutput("stale_note"),
          shiny::h4("Sonuç tablosu"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("table") else shiny::tableOutput("table_fallback"),
          shiny::h4("AUC ısı haritası"),
          shiny::plotOutput("heatmap", height = "320px"),
          shiny::h4("Dayanıklılık (Sensitivity)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("sensitivity_table") else shiny::tableOutput("sensitivity_fallback")
        ),
        shiny::tabPanel(
          "Gen Kapsama",
          shiny::h4("Gen kapsama ve eksik gen detayı"),
          shiny::tags$p(shiny::em("Her imza × kohort için bulunan/toplam gen, % kapsama ve missing_genes.")),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("coverage_table") else shiny::tableOutput("coverage_fallback")
        ),
        shiny::tabPanel(
          "ROC eğrileri",
          shiny::selectInput(
            "roc_dataset", "Veri seti",
            choices = c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
          ),
          shiny::plotOutput("roc_plot", height = "420px")
        ),
        shiny::tabPanel(
          "DeLong",
          shiny::verbatimTextOutput("cmp_status"),
          shiny::tags$p(shiny::em("Tek karşılaştırmada p_adj_bh ham p değerine eşittir.")),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("cmp_table") else shiny::tableOutput("cmp_fallback")
        ),
        shiny::tabPanel(
          "Özet",
          shiny::selectInput(
            "forest_sig", "Forest imzası",
            choices = c("Cytokine-5", "IFNG6", "Random3")
          ),
          shiny::h4("Ağırlıklı Özet (n + ters-varyans, Q, I2)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("summary_table") else shiny::tableOutput("summary_fallback"),
          shiny::h4("AUC forest"),
          shiny::plotOutput("forest_plot", height = "320px"),
          shiny::h4("Jaccard (aktif imzalar)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("jaccard_table") else shiny::tableOutput("jaccard_fallback"),
          shiny::h4("AUC korelasyonu (veri setleri boyunca)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("corr_table") else shiny::tableOutput("corr_fallback")
        ),
        shiny::tabPanel(
          "R kodu",
          shiny::h4("Bu analizin R kodunu kopyala"),
          shiny::tags$p("Son çalıştırılan snapshot parametrelerinden üretilen kod:"),
          shiny::textAreaInput(
            "r_code_text",
            label = NULL,
            value = "",
            rows = 18,
            width = "100%"
          ),
          shiny::tags$button(
            "Kodu kopyala",
            id = "copy_r_code",
            class = "btn btn-default",
            type = "button",
            onclick = "sigbenchrCopyRcode();"
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
        paste0("Hata: ", conditionMessage(ds))
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
        "Eklendi: ", nm,
        " | gen=", ds$info$n_genes,
        " | ornek=", ds$info$n_samples,
        " | norm=", ds$info$normalization,
        " | scale=", ds$info$scale
      )
      if (length(warn_msgs) > 0L) {
        msg <- paste0(msg, "\nUyarilar:\n- ", paste(warn_msgs, collapse = "\n- "))
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
      shiny::validate(shiny::need(nm %in% names(signatures_demo), paste0("Bilinmeyen imza: ", nm)))
      sig <- signatures_demo[[nm]]
      sig$direction <- direction
      sigs[[nm]] <- sig
    }
    custom_txt <- if (is.null(custom_genes)) "" else custom_genes
    if (nzchar(trimws(custom_txt))) {
      genes <- SigBenchR::parse_genes(custom_txt)
      shiny::validate(shiny::need(length(genes) >= 1, "Özel imza için en az bir gen girin."))
      # derived_from belirtilmez -> independent = NA (fail-closed)
      sigs[["Custom"]] <- list(genes = genes, direction = direction)
    }
    shiny::validate(shiny::need(length(sigs) >= 1, "En az bir hazır imza seçin veya özel imza ekleyin."))
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
      paste0("Bilinmeyen veri seti: ", paste(miss_ds, collapse = ", "))
    ))
    ds <- pool[params$datasets]

    shiny::withProgress(message = "Benchmark çalışıyor...", value = 0.2, {
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
      shiny::tags$p(
        shiny::em("Ayarlar değişti; güncel sonuç için 'Analizi Çalıştır'a basın."),
        style = "color:#8a6d3b;"
      )
    } else {
      NULL
    }
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
    suppressWarnings(
      SigBenchR::summarize_auc(snap$res, method = hm)
    )
  })

  shiny::observe({
    snap <- results()
    if (is.null(snap)) {
      shiny::updateTextAreaInput(
        session, "r_code_text",
        value = "# Analizi çalıştırdıktan sonra tekrarlanabilir R kodu burada görünür."
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
        compute_null = p$compute_null
      ),
      error = function(e) paste0("# Kod üretilemedi: ", conditionMessage(e))
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
        stop("HTML rapor için 'rmarkdown' paketini kurun.", call. = FALSE)
      }
      snap <- results()
      SigBenchR::export_report(
        results = snap$res,
        summary_tbl = summary_res(),
        sensitivity_tbl = sens(),
        params_list = list(
          methods = snap$params$methods,
          direction = snap$params$direction,
          coverage_threshold = snap$params$coverage_threshold,
          sensitivity_threshold = snap$params$sensitivity_threshold,
          n_boot = snap$params$n_boot,
          n_perm = snap$params$n_perm,
          seed = snap$params$seed,
          sd_handling = snap$params$sd_handling,
          match_expression = snap$params$match_expression,
          heatmap_method = snap$params$heatmap_method,
          forest_signature = snap$params$forest_signature
        ),
        output_file = file
      )
    }
  )

  output$status <- shiny::renderText({
    if (is.null(input$run) || identical(as.integer(input$run), 0L)) {
      return("Analizi çalıştırmak için düğmeye basın.")
    }
    if (is.null(results())) {
      return("Analizi çalıştırmak için düğmeye basın.")
    }
    df <- results()$res
    paste0(
      "Satır: ", nrow(df),
      " | Düşük kapsama: ", sum(df$low_coverage, na.rm = TRUE),
      " | Küçük kohort: ", sum(df$small_cohort, na.rm = TRUE),
      " | var_zero: ", sum(df$var_zero, na.rm = TRUE)
    )
  })

  .dt_round <- function(df) {
    if (!requireNamespace("DT", quietly = TRUE)) {
      return(df)
    }
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    dt <- DT::datatable(df, options = list(scrollX = TRUE, pageLength = 8))
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
      "AUC yok."
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
      shiny::validate(shiny::need(FALSE, "ROC için önce analizi çalıştırın."))
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
        "DeLong için iki farklı imza seçin (İmza 1 ≠ İmza 2)."
      )
    )
    shiny::withProgress(message = "DeLong karşılaştırması...", value = 0.4, {
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
        shiny::need(length(genes1) >= 1, "İmza 1 genleri bulunamadı."),
        shiny::need(length(genes2) >= 1, "İmza 2 genleri bulunamadı.")
      )

      cmp_pool <- active_dataset_pool()
      shiny::validate(shiny::need(
        input$cmp_dataset %in% names(cmp_pool),
        paste0("Veri seti bulunamadı: ", input$cmp_dataset)
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
        paste0("İmza 1 skorlanamadı: ", if (inherits(s1, "error")) conditionMessage(s1) else "")
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
        paste0("İmza 2 skorlanamadı: ", if (inherits(s2, "error")) conditionMessage(s2) else "")
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
        paste0("DeLong karşılaştırma hatası: ", if (inherits(cmp, "error")) conditionMessage(cmp) else "")
      ))
      cmp
    })
  }, ignoreInit = TRUE)

  output$cmp_status <- shiny::renderText({
    if (identical(input$cmp_sig1, input$cmp_sig2)) {
      return("Uyarı: İmza 1 ve İmza 2 aynı olamaz.")
    }
    if (is.null(input$run_cmp) || identical(as.integer(input$run_cmp), 0L) || is.null(cmp_res())) {
      return("Aynı veri setinde iki imza seçip DeLong karşılaştırın.")
    }
    r <- cmp_res()
    sprintf(
      "p=%.4g | p_adj_bh=%.4g (tek karşılaştırmada ham p ile eşit) | AUC1=%.3f | AUC2=%.3f",
      r$p_value, r$p_adj_bh, r$auc1, r$auc2
    )
  })
  output$cmp_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(cmp_res()); .dt_round(cmp_res()) })
  } else NULL
  output$cmp_fallback <- shiny::renderTable({ shiny::req(cmp_res()); cmp_res() })

  output$summary_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(summary_res()); .dt_round(summary_res()) })
  } else NULL
  output$summary_fallback <- shiny::renderTable({ shiny::req(summary_res()); summary_res() })

  output$forest_plot <- shiny::renderPlot({
    shiny::req(results())
    snap <- results()
    shiny::validate(shiny::need(
      any(snap$res$signature == input$forest_sig &
            snap$res$method == input$heatmap_method &
            !(snap$res$skipped %in% TRUE)),
      "Seçilen imza/yöntem için AUC yok."
    ))
    SigBenchR::plot_auc_forest(
      snap$res,
      signature = input$forest_sig,
      method = input$heatmap_method,
      summary = summary_res()
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
        "Korelasyon için en az iki imza çalıştırın."
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
