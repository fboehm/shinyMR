library(shiny)
library(shinyjs)
library(TwoSampleMR)
library(dplyr)
library(rmarkdown)

# Helper function: Read uploaded file safely
read_uploaded_file <- function(file, sep) {
  req(file)
  tryCatch({
    data <- read.table(file$datapath, header = TRUE, sep = sep, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
    
    if (ncol(data) < 8) {
      stop("Error: The uploaded file does not have enough columns. Please check the file format.")
    }
    return(data)
  }, error = function(e) {
    showNotification("Error reading file. Please check the file format.", type = "error")
    return(NULL)
  })
}

# UI Definition
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Two-Sample Mendelian Randomization (MR) with Sensitivity Analyses"),
  sidebarLayout(
    sidebarPanel(
      h4("Upload GWAS Summary Statistics"),
      fileInput("exposure_file", "Upload Exposure GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      fileInput("outcome_file", "Upload Outcome GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      selectInput("sep", "File Separator:", choices = c("Tab" = "\t", "Comma" = ","), selected = "\t"),
      numericInput("pval_threshold", "P-value threshold for SNP selection:", value = 5e-8, step = 1e-8),
      actionButton("update_threshold", "Update P-Value Threshold"),
      hr(),
      h4("Column Mapping for Exposure"),
      uiOutput("exposure_columns"),
      h4("Column Mapping for Outcome"),
      uiOutput("outcome_columns"),
      hr(),
      actionButton("run_analysis", "Run MR & Sensitivity Analyses"),
      actionButton("reset_uploads", "Reset Uploads", style = "background-color: red; color: white;"),
      downloadButton("download_report", "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Harmonized Data", tableOutput("harmonized_data")),
        tabPanel("MR Results", tableOutput("mr_results")),
        tabPanel("Heterogeneity Test", tableOutput("heterogeneity_results")),
        tabPanel("Pleiotropy Test", tableOutput("pleiotropy_results")),
        tabPanel("Single SNP Analysis", tableOutput("single_snp_results")),
        tabPanel("Leave-One-Out Analysis", tableOutput("leave_one_out_results")),
        tabPanel("Plots", 
                 plotOutput("mr_plot"), 
                 plotOutput("funnel_plot"), 
                 plotOutput("loo_forest_plot"))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(exposure_data = NULL, outcome_data = NULL, pval_threshold = 5e-8)

  observeEvent(input$exposure_file, {
    rv$exposure_data <- read_uploaded_file(input$exposure_file, input$sep)
  })

  observeEvent(input$outcome_file, {
    rv$outcome_data <- read_uploaded_file(input$outcome_file, input$sep)
  })

  observeEvent(input$update_threshold, {
    rv$pval_threshold <- input$pval_threshold
  })

  harmonized_data <- eventReactive(input$run_analysis, {
    req(rv$exposure_data, rv$outcome_data)
    exposure_filtered <- rv$exposure_data %>% filter(get(input$exposure_pval) < rv$pval_threshold)
    harmonise_data(exposure_dat = exposure_filtered, outcome_dat = rv$outcome_data)
  })

  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr(harmonized_data())
  })

  heterogeneity_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr_heterogeneity(harmonized_data())
  })

  pleiotropy_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr_pleiotropy_test(harmonized_data())
  })

  single_snp_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr_singlesnp(harmonized_data())
  })

  leave_one_out_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr_leaveoneout(harmonized_data())
  })

  output$harmonized_data <- renderTable({ req(harmonized_data()); head(harmonized_data()) })
  output$mr_results <- renderTable({ req(mr_results()); mr_results() })
  output$heterogeneity_results <- renderTable({ req(heterogeneity_results()); heterogeneity_results() })
  output$pleiotropy_results <- renderTable({ req(pleiotropy_results()); pleiotropy_results() })
  output$single_snp_results <- renderTable({ req(single_snp_results()); single_snp_results() })
  output$leave_one_out_results <- renderTable({ req(leave_one_out_results()); leave_one_out_results() })

  output$mr_plot <- renderPlot({ req(mr_results(), harmonized_data()); mr_scatter_plot(mr_results(), harmonized_data()) })
  output$funnel_plot <- renderPlot({ req(single_snp_results()); mr_funnel_plot(single_snp_results()) })
  output$loo_forest_plot <- renderPlot({ req(leave_one_out_results()); mr_leaveoneout_plot(leave_one_out_results()) })

  observeEvent(input$reset_uploads, {
    rv$exposure_data <- NULL
    rv$outcome_data <- NULL
    showNotification("Uploads reset successfully!", type = "message")
  })
}

shinyApp(ui = ui, server = server)
