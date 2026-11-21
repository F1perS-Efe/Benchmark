# SigBenchR Shiny MVP (Faz 1) - tek dosya

options(shiny.maxRequestSize = 25 * 1024^2) # 25 MB

ui <- shiny::fluidPage(
  shiny::titlePanel("SigBenchR - Imza Dogrulama (Faz 1)"),
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
      shiny::radioButtons(
        "method",
        "Skorlama yontemi",
        choices = c("mean", "zscore"),
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
        "AUC her veri seti icinde ayri hesaplanir."
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
      shiny::plotOutput("heatmap", height = "360px")
    )
  )
)

server <- function(input, output, session) {
  data(datasets_demo, package = "SigBenchR", envir = environment())
  data(signatures_demo, package = "SigBenchR", envir = environment())

  results <- shiny::eventReactive(input$run, {
    shiny::req(length(input$datasets) > 0)

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
      methods = input$method,
      direction = input$direction
    )
  })

  output$status <- shiny::renderText({
    if (is.null(results())) {
      return("Analizi calistirmak icin dugmeye basin.")
    }
    df <- results()
    paste0(
      "Satir: ", nrow(df),
      " | Atlanan: ", sum(df$skipped, na.rm = TRUE),
      " | Dusuk kapsama: ", sum(df$low_coverage, na.rm = TRUE)
    )
  })

  output$table <- if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT({
      shiny::req(results())
      df <- results()
      cols <- c(
        "dataset", "signature", "method", "n", "n_pos", "n_neg",
        "coverage_pct", "auc", "normalization", "scale",
        "cancer", "treatment", "low_coverage", "skipped", "skip_reason"
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
    cols <- c("dataset", "signature", "method", "n", "coverage_pct", "auc", "skipped")
    df[, intersect(cols, names(df)), drop = FALSE]
  })

  output$heatmap <- shiny::renderPlot({
    shiny::req(results())
    df <- results()
    shiny::validate(shiny::need(any(!df$skipped & is.finite(df$auc)), "Cizilecek AUC yok."))
    SigBenchR::plot_auc_heatmap(df, method = input$method)
  })
}

shiny::shinyApp(ui, server)
