# SigBenchR Shiny (Faz 3A) - tek dosya

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
  shiny::titlePanel("SigBenchR - Imza Dogrulama (Faz 3A)"),
  shiny::tags$p(
    shiny::em(
      "Arastirma amaclidir; klinik karar araci degildir. ",
      "Kimliklendirilebilir hasta verisi yuklemeyin."
    )
  ),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::selectInput(
        "sig_choice", "Hazir imza",
        choices = c("Cytokine-5", "IFNG6", "Random3", "Ozel (yapistir)" = "custom"),
        selected = "Cytokine-5"
      ),
      shiny::conditionalPanel(
        "input.sig_choice == 'custom'",
        shiny::textAreaInput("custom_genes", "Gen listesi", rows = 3)
      ),
      shiny::checkboxGroupInput(
        "methods", "Skorlama yontemleri",
        choices = c("mean", "zscore", "median", "rank"),
        selected = c("mean", "zscore")
      ),
      shiny::radioButtons(
        "heatmap_method", "Isi haritasi / ROC yontemi",
        choices = c("mean", "zscore", "median", "rank"), selected = "mean"
      ),
      shiny::radioButtons(
        "direction", "Imza yonu",
        choices = c("Yuksek=yanitli (<)" = "<", "Dusuk=yanitli (>)" = ">"),
        selected = "<"
      ),
      shiny::sliderInput("coverage_threshold", "Kapsama esigi", 0.5, 1, 0.70, 0.05),
      shiny::sliderInput("sensitivity_threshold", "Dayaniklilik esigi (keyfi)", 0.05, 0.3, 0.10, 0.01),
      shiny::checkboxInput("compute_pr", "PR-AUC + bootstrap GA", TRUE),
      shiny::numericInput("n_boot", "PR bootstrap (n_boot)", value = 200, min = 0, max = 2000, step = 50),
      shiny::checkboxInput("compute_null", "Null permutasyon", FALSE),
      shiny::numericInput("n_perm", "Null B", value = 50, min = 10, max = 500, step = 10),
      shiny::numericInput("seed", "Seed", value = 42, min = 1, step = 1),
      shiny::checkboxGroupInput(
        "datasets", "Veri setleri",
        choices = c("demo_A", "demo_B"), selected = c("demo_A", "demo_B")
      ),
      shiny::actionButton("run", "Analizi calistir", class = "btn-primary"),
      shiny::downloadButton("download_report", "Tum sonuclari HTML rapor olarak indir"),
      shiny::hr(),
      shiny::h5("DeLong: iki imza (ayni veri seti)"),
      shiny::selectInput("cmp_dataset", "Veri seti", choices = c("demo_A", "demo_B")),
      shiny::selectInput("cmp_sig1", "Imza 1", choices = c("Cytokine-5", "IFNG6", "Random3")),
      shiny::selectInput("cmp_sig2", "Imza 2", choices = c("IFNG6", "Cytokine-5", "Random3")),
      shiny::actionButton("run_cmp", "DeLong karsilastir")
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Sonuclar",
          shiny::verbatimTextOutput("status"),
          shiny::h4("Sonuc tablosu"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("table") else shiny::tableOutput("table_fallback"),
          shiny::h4("AUC isi haritasi"),
          shiny::plotOutput("heatmap", height = "320px"),
          shiny::h4("Dayaniklilik (Sensitivity)"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("sensitivity_table") else shiny::tableOutput("sensitivity_fallback")
        ),
        shiny::tabPanel(
          "ROC egirileri",
          shiny::selectInput("roc_dataset", "Veri seti", choices = c("demo_A", "demo_B")),
          shiny::plotOutput("roc_plot", height = "420px")
        ),
        shiny::tabPanel(
          "DeLong",
          shiny::verbatimTextOutput("cmp_status"),
          if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("cmp_table") else shiny::tableOutput("cmp_fallback")
        ),
        shiny::tabPanel(
          "Ozet",
          shiny::selectInput(
            "forest_sig", "Forest imzasi",
            choices = c("Cytokine-5", "IFNG6", "Random3")
          ),
          shiny::h4("Agirlikli AUC ozeti (n + ters-varyans, Q, I2)"),
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
          shiny::tags$p("Secili parametrelerden uretilen, dogrudan calistirilabilir kod:"),
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

  build_sigs <- function() {
    if (identical(input$sig_choice, "custom")) {
      genes <- SigBenchR::parse_genes(input$custom_genes)
      shiny::validate(shiny::need(length(genes) >= 1, "Ozel imza icin gen girin."))
      list(Custom = list(genes = genes, direction = input$direction))
    } else {
      sig <- signatures_demo[[input$sig_choice]]
      sig$direction <- input$direction
      setNames(list(sig), input$sig_choice)
    }
  }

  results <- shiny::eventReactive(input$run, {
    shiny::req(length(input$datasets) > 0, length(input$methods) > 0)
    sigs <- build_sigs()
    ds <- datasets_demo[input$datasets]

    shiny::withProgress(message = "Benchmark calisiyor...", value = 0.2, {
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
    sig_name <- if (identical(input$sig_choice, "custom")) NULL else input$sig_choice
    genes <- if (identical(input$sig_choice, "custom")) {
      SigBenchR::parse_genes(input$custom_genes)
    } else {
      NULL
    }
    code <- tryCatch(
      SigBenchR::generate_r_code(
        signature_name = sig_name,
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
      error = function(e) paste0("# Kod uretilemedi: ", conditionMessage(e))
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
        stop("HTML rapor icin 'rmarkdown' paketini kurun.", call. = FALSE)
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
    if (is.null(results())) return("Analizi calistirmak icin dugmeye basin.")
    df <- results()
    paste0(
      "Satir: ", nrow(df),
      " | Dusuk kapsama: ", sum(df$low_coverage, na.rm = TRUE),
      " | Kucuk kohort: ", sum(df$small_cohort, na.rm = TRUE),
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
    shiny::validate(shiny::need(any(df$method == input$heatmap_method & !df$skipped), "AUC yok."))
    SigBenchR::plot_auc_heatmap(df, method = input$heatmap_method, coverage_threshold = input$coverage_threshold)
  })

  output$sensitivity_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({ shiny::req(sens()); DT::datatable(sens()) })
  } else NULL
  output$sensitivity_fallback <- shiny::renderTable({ shiny::req(sens()); sens() })

  output$roc_plot <- shiny::renderPlot({
    shiny::req(input$roc_dataset)
    # Benchmark'taki imzalar veya hazir demolar
    sigs <- if (!is.null(results())) {
      # Son calistirmadaki imza adlari
      sn <- unique(results()$signature)
      if (all(sn %in% names(signatures_demo))) signatures_demo[sn] else build_sigs()
    } else signatures_demo
    SigBenchR::plot_roc_curves(
      datasets = datasets_demo,
      signatures = sigs,
      dataset_name = input$roc_dataset,
      method = input$heatmap_method,
      direction = input$direction
    )
  })

  cmp_res <- shiny::eventReactive(input$run_cmp, {
    shiny::withProgress(message = "DeLong karsilastirmasi...", value = 0.4, {
      s1 <- SigBenchR::score_for_dataset(
        datasets_demo[[input$cmp_dataset]],
        genes = signatures_demo[[input$cmp_sig1]]$genes,
        method = input$heatmap_method,
        direction = input$direction
      )
      s2 <- SigBenchR::score_for_dataset(
        datasets_demo[[input$cmp_dataset]],
        genes = signatures_demo[[input$cmp_sig2]]$genes,
        method = input$heatmap_method,
        direction = input$direction
      )
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
      )
    })
  })

  output$cmp_status <- shiny::renderText({
    if (is.null(cmp_res())) return("Ayni veri setinde iki imza secip DeLong karsilastir.")
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
            !results()$skipped),
      "Secilen imza/yontem icin AUC yok."
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
        length(unique(results()$signature[!results()$skipped])) >= 2,
        "Korelasyon icin en az iki imza calistirin (ornek: hazir imzalari tek tek degil, coklu benchmark)."
      ))
      # Tek imza seciliyse demo imzalarinin hepsiyle yeniden degil; mevcut sonuclardan
      DT::datatable(SigBenchR::auc_correlation(results(), method = input$heatmap_method))
    })
  } else NULL
  output$corr_fallback <- shiny::renderTable({
    shiny::req(results())
    SigBenchR::auc_correlation(results(), method = input$heatmap_method)
  })
}

shiny::shinyApp(ui, server)
