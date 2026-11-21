# SigBenchR Shiny (Faz 2B) - tek dosya

options(shiny.maxRequestSize = 25 * 1024^2)

ui <- shiny::fluidPage(
  shiny::titlePanel("SigBenchR - Imza Dogrulama (Faz 2B)"),
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
      shiny::checkboxInput("compute_null", "Null permütasyon", FALSE),
      shiny::numericInput("n_perm", "Null B", value = 50, min = 10, max = 500, step = 10),
      shiny::numericInput("seed", "Seed", value = 42, min = 1, step = 1),
      shiny::checkboxGroupInput(
        "datasets", "Veri setleri",
        choices = c("demo_A", "demo_B"), selected = c("demo_A", "demo_B")
      ),
      shiny::actionButton("run", "Analizi calistir", class = "btn-primary"),
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
        "pr_auc", "pr_baseline", "pr_ci_low", "pr_ci_high", "null_p",
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
}

shiny::shinyApp(ui, server)
