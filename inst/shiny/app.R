# SigBenchR Shiny (Faz 2A) - tek dosya

options(shiny.maxRequestSize = 25 * 1024^2) # 25 MB

ui <- shiny::fluidPage(
  shiny::titlePanel("SigBenchR - Imza Dogrulama (Faz 2A)"),
  shiny::tags$p(
    shiny::em(
      "Arastirma amaclidir; klinik karar araci degildir. ",
      "Kimliklendirilebilir hasta verisi yuklemeyin. ",
      "Yuklenen veriler diske kaydedilmez."
    )
  ),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::selectInput(
        "sig_choice",
        "Hazir imza",
        choices = c(
          "Cytokine-5" = "Cytokine-5",
          "IFNG6" = "IFNG6",
          "Random3" = "Random3",
          "Ozel (yapistir)" = "custom"
        ),
        selected = "Cytokine-5"
      ),
      shiny::conditionalPanel(
        condition = "input.sig_choice == 'custom'",
        shiny::textAreaInput(
          "custom_genes",
          "Gen listesi (virgul / satir)",
          rows = 4,
          placeholder = "IFNG, CXCL9, CXCL13"
        )
      ),
      shiny::checkboxGroupInput(
        "methods",
        "Skorlama yontemleri",
        choices = c("mean", "zscore", "median", "rank"),
        selected = c("mean", "zscore", "median", "rank")
      ),
      shiny::radioButtons(
        "heatmap_method",
        "Isi haritasi yontemi",
        choices = c("mean", "zscore", "median", "rank"),
        selected = "mean"
      ),
      shiny::radioButtons(
        "direction",
        "Imza yonu",
        choices = c(
          "Yuksek skor = yanitli (<)" = "<",
          "Dusuk skor = yanitli (>)" = ">"
        ),
        selected = "<"
      ),
      shiny::sliderInput(
        "coverage_threshold",
        "Kapsama esigi",
        min = 0.5, max = 1, value = 0.70, step = 0.05
      ),
      shiny::sliderInput(
        "sensitivity_threshold",
        "Dayaniklilik AUC-aralik esigi (keyfi)",
        min = 0.05, max = 0.30, value = 0.10, step = 0.01
      ),
      shiny::tags$small(
        "0.10 esigi keyfidir; bilimsel bir standart degildir. ",
        "Sonuclara bakildiktan sonra degistirilmemelidir."
      ),
      shiny::checkboxGroupInput(
        "datasets",
        "Veri setleri",
        choices = c("demo_A", "demo_B"),
        selected = c("demo_A", "demo_B")
      ),
      shiny::actionButton("run", "Analizi calistir", class = "btn-primary"),
      shiny::hr(),
      shiny::tags$small(
        "Farkli normalizasyonlu veri setlerinin ham ifadeleri birlestirilmez. ",
        "AUC her veri seti icinde ayri hesaplanir. ",
        "Gri isi-haritasi hucreleri dusuk kapsama uyarisidir."
      )
    ),
    shiny::mainPanel(
      shiny::verbatimTextOutput("status"),
      shiny::h4("Sonuc tablosu"),
      if (requireNamespace("DT", quietly = TRUE)) {
        DT::DTOutput("table")
      } else {
        shiny::tableOutput("table_fallback")
      },
      shiny::h4("AUC isi haritasi"),
      shiny::plotOutput("heatmap", height = "360px"),
      shiny::h4("Dayaniklilik (Sensitivity)"),
      shiny::tags$p(
        shiny::em(
          "fragile = TRUE: yontemler arasi AUC araligi esigi asar. ",
          "Esik keyfidir (varsayilan 0.10)."
        )
      ),
      if (requireNamespace("DT", quietly = TRUE)) {
        DT::DTOutput("sensitivity_table")
      } else {
        shiny::tableOutput("sensitivity_fallback")
      }
    )
  )
)

server <- function(input, output, session) {
  data(datasets_demo, package = "SigBenchR", envir = environment())
  data(signatures_demo, package = "SigBenchR", envir = environment())

  results <- shiny::eventReactive(input$run, {
    shiny::req(length(input$datasets) > 0)
    shiny::validate(shiny::need(length(input$methods) >= 1, "En az bir skorlama yontemi secin."))

    if (identical(input$sig_choice, "custom")) {
      genes <- SigBenchR::parse_genes(input$custom_genes)
      shiny::validate(shiny::need(length(genes) >= 1, "Ozel imza icin en az bir gen girin."))
      sigs <- list(Custom = list(genes = genes, direction = input$direction))
    } else {
      sig <- signatures_demo[[input$sig_choice]]
      sig$direction <- input$direction
      sigs <- list()
      sigs[[input$sig_choice]] <- sig
    }

    ds <- datasets_demo[input$datasets]
    SigBenchR::run_benchmark(
      signatures = sigs,
      datasets = ds,
      methods = input$methods,
      direction = input$direction,
      coverage_threshold = input$coverage_threshold
    )
  })

  sens <- shiny::reactive({
    shiny::req(results())
    SigBenchR::check_sensitivity(
      results(),
      threshold = input$sensitivity_threshold
    )
  })

  output$status <- shiny::renderText({
    if (is.null(results())) {
      return("Analizi calistirmak icin dugmeye basin.")
    }
    df <- results()
    s <- sens()
    paste0(
      "Satir: ", nrow(df),
      " | Atlanan: ", sum(df$skipped, na.rm = TRUE),
      " | Dusuk kapsama: ", sum(df$low_coverage, na.rm = TRUE),
      " | Kirilgan (fragile): ", sum(s$fragile, na.rm = TRUE)
    )
  })

  output$table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      df <- results()
      cols <- c(
        "dataset", "signature", "method", "n", "n_pos", "n_neg",
        "coverage_label", "missing_genes", "low_coverage", "auc",
        "normalization", "scale", "cancer", "treatment",
        "skipped", "skip_reason"
      )
      cols <- intersect(cols, names(df))
      DT::datatable(df[, cols, drop = FALSE], options = list(pageLength = 10))
    })
  } else {
    NULL
  }

  output$table_fallback <- shiny::renderTable({
    shiny::req(results())
    df <- results()
    cols <- c(
      "dataset", "signature", "method", "n", "coverage_label",
      "missing_genes", "low_coverage", "auc", "skipped"
    )
    df[, intersect(cols, names(df)), drop = FALSE]
  })

  output$heatmap <- shiny::renderPlot({
    shiny::req(results())
    df <- results()
    shiny::validate(shiny::need(
      any(df$method == input$heatmap_method & !df$skipped & is.finite(df$auc)),
      "Secilen isi-haritasi yontemi icin cizilecek AUC yok."
    ))
    SigBenchR::plot_auc_heatmap(
      df,
      method = input$heatmap_method,
      coverage_threshold = input$coverage_threshold
    )
  })

  output$sensitivity_table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(sens())
      DT::datatable(sens(), options = list(pageLength = 10))
    })
  } else {
    NULL
  }

  output$sensitivity_fallback <- shiny::renderTable({
    shiny::req(sens())
    sens()
  })
}

shiny::shinyApp(ui, server)
