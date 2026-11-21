# SigBenchR Shiny - tek dosya

options(shiny.maxRequestSize = 25 * 1024^2)

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
      "Araştırma amaçlıdır; klinik karar aracı değildir. ",
      "Kimliklendirilebilir hasta verisi yüklemeyin."
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
        "heatmap_method", "Isı Haritası / ROC yöntemi",
        choices = c("mean", "zscore"), selected = "mean"
      ),
      shiny::radioButtons(
        "direction", "İmza yönü",
        choices = c("Yüksek=yanıtlı (<)" = "<", "Düşük=yanıtlı (>)" = ">"),
        selected = "<"
      ),
      shiny::sliderInput("coverage_threshold", "Kapsama eşiği", 0.5, 1, 0.70, 0.05),
      shiny::sliderInput("sensitivity_threshold", "Dayanıklılık eşiği (keyfi)", 0.05, 0.3, 0.10, 0.01),
      shiny::checkboxInput("compute_pr", "PR-AUC + bootstrap GA", TRUE),
      shiny::numericInput("n_boot", "PR bootstrap (n_boot)", value = 200, min = 0, max = 2000, step = 50),
      shiny::checkboxInput("compute_null", "Null permütasyon", FALSE),
      shiny::numericInput("n_perm", "Null B", value = 200, min = 10, max = 2000, step = 10),
      shiny::numericInput("seed", "Seed", value = 42, min = 1, step = 1),
      shiny::checkboxGroupInput(
        "datasets", "Veri setleri",
        choices = c("demo_A", "demo_B"), selected = c("demo_A", "demo_B")
      ),
      shiny::actionButton("run", "Analizi çalıştır", class = "btn-primary"),
      shiny::downloadButton("download_report", "Tüm sonuçları HTML rapor olarak indir"),
      shiny::hr(),
      shiny::h5("DeLong: iki imza (aynı veri seti)"),
      shiny::selectInput("cmp_dataset", "Veri seti", choices = c("demo_A", "demo_B")),
      shiny::selectInput("cmp_sig1", "İmza 1", choices = c("Cytokine-5", "IFNG6", "Random3")),
      shiny::selectInput("cmp_sig2", "İmza 2", choices = c("IFNG6", "Cytokine-5", "Random3")),
      shiny::actionButton("run_cmp", "DeLong karşılaştır")
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Sonuçlar",
          shiny::verbatimTextOutput("status"),
          shiny::h4("Sonuç tablosu"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("table") else shiny::tableOutput("table_fallback"),
          shiny::h4("AUC ısı haritası"),
          shiny::plotOutput("heatmap", height = "320px"),
          shiny::h4("Dayanıklılık (Sensitivity)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("sensitivity_table") else shiny::tableOutput("sensitivity_fallback")
        ),
        shiny::tabPanel(
          "ROC eğrileri",
          shiny::selectInput("roc_dataset", "Veri seti", choices = c("demo_A", "demo_B")),
          shiny::plotOutput("roc_plot", height = "420px")
        ),
        shiny::tabPanel(
          "DeLong",
          shiny::verbatimTextOutput("cmp_status"),
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
          shiny::h4("Jaccard (demo imzalar)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("jaccard_table") else shiny::tableOutput("jaccard_fallback"),
          shiny::h4("AUC korelasyonu (veri setleri boyunca)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("corr_table") else shiny::tableOutput("corr_fallback")
        ),
        shiny::tabPanel(
          "R kodu",
          shiny::h4("Bu analizin R kodunu kopyala"),
          shiny::tags$p("Seçili parametrelerden üretilen, doğrudan çalıştırılabilir kod:"),
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
  data(datasets_demo, package = "SigBenchR", envir = environment())
  data(signatures_demo, package = "SigBenchR", envir = environment())

  active_sig_names <- shiny::reactive({
    chosen <- as.character(input$sig_choice)
    chosen <- chosen[nzchar(chosen)]
    custom_txt <- if (is.null(input$custom_genes)) "" else input$custom_genes
    if (nzchar(trimws(custom_txt))) {
      chosen <- c(chosen, "Custom")
    }
    unique(chosen)
  })

  build_sigs <- function() {
    chosen <- as.character(input$sig_choice)
    chosen <- chosen[nzchar(chosen)]
    sigs <- list()
    for (nm in chosen) {
      shiny::validate(shiny::need(nm %in% names(signatures_demo), paste0("Bilinmeyen imza: ", nm)))
      sig <- signatures_demo[[nm]]
      sig$direction <- input$direction
      sigs[[nm]] <- sig
    }
    custom_txt <- if (is.null(input$custom_genes)) "" else input$custom_genes
    if (nzchar(trimws(custom_txt))) {
      genes <- SigBenchR::parse_genes(custom_txt)
      shiny::validate(shiny::need(length(genes) >= 1, "Özel imza için en az bir gen girin."))
      sigs[["Custom"]] <- list(genes = genes, direction = input$direction)
    }
    shiny::validate(shiny::need(length(sigs) >= 1, "En az bir hazır imza seçin veya özel imza ekleyin."))
    sigs
  }

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

  results <- shiny::eventReactive(input$run, {
    shiny::req(length(input$datasets) > 0, length(input$methods) > 0)
    sigs <- build_sigs()
    ds <- datasets_demo[input$datasets]

    shiny::withProgress(message = "Benchmark çalışıyor...", value = 0.2, {
      shiny::incProgress(0.3, detail = "ROC / PR / null")
      out <- SigBenchR::run_benchmark(
        signatures = sigs,
        datasets = ds,
        methods = input$methods,
        direction = input$direction,
        coverage_threshold = input$coverage_threshold,
        compute_pr = isTRUE(input$compute_pr),
        n_boot = as.integer(input$n_boot),
        compute_null = isTRUE(input$compute_null),
        n_perm = as.integer(input$n_perm),
        seed = as.integer(input$seed)
      )
      shiny::incProgress(0.5)
      out
    })
  })

  sens <- shiny::reactive({
    shiny::req(results())
    SigBenchR::check_sensitivity(results(), threshold = input$sensitivity_threshold)
  })

  summary_res <- shiny::reactive({
    shiny::req(results())
    suppressWarnings(
      SigBenchR::summarize_auc(results(), method = input$heatmap_method)
    )
  })

  shiny::observe({
    demo_names <- as.character(input$sig_choice)
    demo_names <- demo_names[nzchar(demo_names)]
    custom_txt <- if (is.null(input$custom_genes)) "" else input$custom_genes
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
        methods = input$methods,
        direction = input$direction,
        coverage_threshold = input$coverage_threshold,
        n_boot = as.integer(input$n_boot),
        n_perm = as.integer(input$n_perm),
        seed = as.integer(input$seed),
        dataset_names = input$datasets,
        compute_pr = isTRUE(input$compute_pr),
        compute_null = isTRUE(input$compute_null)
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
      SigBenchR::export_report(
        results = results(),
        summary_tbl = summary_res(),
        sensitivity_tbl = sens(),
        params_list = list(
          methods = input$methods,
          direction = input$direction,
          coverage_threshold = input$coverage_threshold,
          sensitivity_threshold = input$sensitivity_threshold,
          n_boot = as.integer(input$n_boot),
          n_perm = as.integer(input$n_perm),
          seed = as.integer(input$seed),
          heatmap_method = input$heatmap_method,
          forest_signature = input$forest_sig
        ),
        output_file = file
      )
    }
  )

  output$status <- shiny::renderText({
    if (is.null(results())) return("Analizi çalıştırmak için düğmeye basın.")
    df <- results()
    paste0(
      "Satır: ", nrow(df),
      " | Düşük kapsama: ", sum(df$low_coverage, na.rm = TRUE),
      " | Küçük kohort: ", sum(df$small_cohort, na.rm = TRUE),
      " | var_zero: ", sum(df$var_zero, na.rm = TRUE)
    )
  })

  output$table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      df <- results()
      cols <- c(
        "dataset", "signature", "method", "n", "n_pos", "n_neg",
        "coverage_label", "missing_genes", "low_coverage",
        "auc", "auc_ci_low", "auc_ci_high", "auc_var", "var_zero", "small_cohort",
        "pr_auc", "pr_baseline", "pr_ci_low", "pr_ci_high",
        "null_p", "null_p_adj_bh", "independent",
        "skipped", "skip_reason"
      )
      DT::datatable(df[, intersect(cols, names(df)), drop = FALSE], options = list(scrollX = TRUE, pageLength = 8))
    })
  } else NULL

  output$table_fallback <- shiny::renderTable({
    shiny::req(results())
    results()[, intersect(c("dataset", "signature", "method", "auc", "pr_auc", "null_p", "low_coverage"), names(results())), drop = FALSE]
  })

  output$heatmap <- shiny::renderPlot({
    shiny::req(results())
    df <- results()
    shiny::validate(shiny::need(
      any(df$method == input$heatmap_method & !(df$skipped %in% TRUE)),
      "AUC yok."
    ))
    SigBenchR::plot_auc_heatmap(df, method = input$heatmap_method, coverage_threshold = input$coverage_threshold)
  })

  output$sensitivity_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(sens()); DT::datatable(sens()) })
  } else NULL
  output$sensitivity_fallback <- shiny::renderTable({ shiny::req(sens()); sens() })

  output$roc_plot <- shiny::renderPlot({
    shiny::req(input$roc_dataset)
    sigs <- if (!is.null(results())) {
      sn <- unique(results()$signature)
      if (all(sn %in% names(signatures_demo))) {
        signatures_demo[sn]
      } else {
        build_sigs()
      }
    } else {
      signatures_demo
    }
    SigBenchR::plot_roc_curves(
      datasets = datasets_demo,
      signatures = sigs,
      dataset_name = input$roc_dataset,
      method = input$heatmap_method,
      direction = input$direction
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

      s1 <- tryCatch(
        SigBenchR::score_for_dataset(
          datasets_demo[[input$cmp_dataset]],
          genes = genes1,
          method = input$heatmap_method,
          direction = input$direction
        ),
        error = function(e) e
      )
      shiny::validate(shiny::need(
        !inherits(s1, "error"),
        paste0("İmza 1 skorlanamadı: ", if (inherits(s1, "error")) conditionMessage(s1) else "")
      ))

      s2 <- tryCatch(
        SigBenchR::score_for_dataset(
          datasets_demo[[input$cmp_dataset]],
          genes = genes2,
          method = input$heatmap_method,
          direction = input$direction
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
  })

  output$cmp_status <- shiny::renderText({
    if (identical(input$cmp_sig1, input$cmp_sig2)) {
      return("Uyarı: İmza 1 ve İmza 2 aynı olamaz.")
    }
    if (is.null(cmp_res())) return("Aynı veri setinde iki imza seçip DeLong karşılaştırın.")
    r <- cmp_res()
    sprintf("p=%.4g | p_adj_bh=%.4g | AUC1=%.3f | AUC2=%.3f", r$p_value, r$p_adj_bh, r$auc1, r$auc2)
  })
  output$cmp_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(cmp_res()); DT::datatable(cmp_res()) })
  } else NULL
  output$cmp_fallback <- shiny::renderTable({ shiny::req(cmp_res()); cmp_res() })

  output$summary_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(summary_res()); DT::datatable(summary_res()) })
  } else NULL
  output$summary_fallback <- shiny::renderTable({ shiny::req(summary_res()); summary_res() })

  output$forest_plot <- shiny::renderPlot({
    shiny::req(results())
    shiny::validate(shiny::need(
      any(results()$signature == input$forest_sig &
            results()$method == input$heatmap_method &
            !(results()$skipped %in% TRUE)),
      "Seçilen imza/yöntem için AUC yok."
    ))
    SigBenchR::plot_auc_forest(
      results(),
      signature = input$forest_sig,
      method = input$heatmap_method,
      summary = summary_res()
    )
  })

  output$jaccard_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      DT::datatable(SigBenchR::jaccard_signatures(signatures_demo))
    })
  } else NULL
  output$jaccard_fallback <- shiny::renderTable({
    SigBenchR::jaccard_signatures(signatures_demo)
  })

  output$corr_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      shiny::validate(shiny::need(
        length(unique(results()$signature[!(results()$skipped %in% TRUE)])) >= 2,
        "Korelasyon için en az iki imza çalıştırın."
      ))
      DT::datatable(SigBenchR::auc_correlation(results(), method = input$heatmap_method))
    })
  } else NULL
  output$corr_fallback <- shiny::renderTable({
    shiny::req(results())
    SigBenchR::auc_correlation(results(), method = input$heatmap_method)
  })
}

shiny::shinyApp(ui, server)
