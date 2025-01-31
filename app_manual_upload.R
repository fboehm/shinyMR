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
    
    # Ensure all rows have the correct number of columns
    if (ncol(data) < 8) {  # Adjust this based on the expected number of columns
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

# Server Logic
server <- function(input, output, session) {
  rv <- reactiveValues(exposure_data = NULL, outcome_data = NULL)
  
  # Read uploaded files
  observeEvent(input$exposure_file, {
    rv$exposure_data <- read_uploaded_file(input$exposure_file, input$sep)
  })
  
  observeEvent(input$outcome_file, {
    rv$outcome_data <- read_uploaded_file(input$outcome_file, input$sep)
  })
  
  # Generate dropdowns for column mapping
  output$exposure_columns <- renderUI({
    req(rv$exposure_data)
    column_names <- names(rv$exposure_data)
    tagList(
      selectInput("exposure_snp", "SNP Column:", choices = column_names),
      selectInput("exposure_beta", "Beta Column:", choices = column_names),
      selectInput("exposure_se", "SE Column:", choices = column_names),
      selectInput("exposure_effect_allele", "Effect Allele Column:", choices = column_names),
      selectInput("exposure_other_allele", "Other Allele Column:", choices = column_names),
      selectInput("exposure_eaf", "EAF Column:", choices = column_names),
      selectInput("exposure_pval", "P-Value Column:", choices = column_names),
      selectInput("exposure_samplesize", "Sample Size Column:", choices = column_names)
    )
  })
  
  output$outcome_columns <- renderUI({
    req(rv$outcome_data)
    column_names <- names(rv$outcome_data)
    tagList(
      selectInput("outcome_snp", "SNP Column:", choices = column_names),
      selectInput("outcome_beta", "Beta Column:", choices = column_names),
      selectInput("outcome_se", "SE Column:", choices = column_names),
      selectInput("outcome_effect_allele", "Effect Allele Column:", choices = column_names),
      selectInput("outcome_other_allele", "Other Allele Column:", choices = column_names),
      selectInput("outcome_eaf", "EAF Column:", choices = column_names),
      selectInput("outcome_pval", "P-Value Column:", choices = column_names),
      selectInput("outcome_samplesize", "Sample Size Column:", choices = column_names)
    )
  })
  
  # Harmonize data
  harmonized_data <- eventReactive(input$run_analysis, {
    req(rv$exposure_data, rv$outcome_data)
    harmonise_data(exposure_dat = rv$exposure_data, outcome_dat = rv$outcome_data)
  })
  
  # MR Analysis
  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr(harmonized_data())
  })
  
  # Sensitivity Analyses
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
  
  # Render Tables
  output$harmonized_data <- renderTable({ req(harmonized_data()); head(harmonized_data()) })
  output$mr_results <- renderTable({ req(mr_results()); mr_results() })
  output$heterogeneity_results <- renderTable({ req(heterogeneity_results()); heterogeneity_results() })
  output$pleiotropy_results <- renderTable({ req(pleiotropy_results()); pleiotropy_results() })
  output$single_snp_results <- renderTable({ req(single_snp_results()); single_snp_results() })
  output$leave_one_out_results <- renderTable({ req(leave_one_out_results()); leave_one_out_results() })
  
  # Generate Plots
  output$mr_plot <- renderPlot({ req(mr_results(), harmonized_data()); mr_scatter_plot(mr_results(), harmonized_data()) })
  output$funnel_plot <- renderPlot({ req(single_snp_results()); mr_funnel_plot(single_snp_results()) })
  output$loo_forest_plot <- renderPlot({ req(leave_one_out_results()); mr_leaveoneout_plot(leave_one_out_results()) })
  
  # Reset uploaded files
  observeEvent(input$reset_uploads, {
    rv$exposure_data <- NULL
    rv$outcome_data <- NULL
    output$harmonized_data <- renderTable(NULL)
    output$mr_results <- renderTable(NULL)
    output$heterogeneity_results <- renderTable(NULL)
    output$pleiotropy_results <- renderTable(NULL)
    output$single_snp_results <- renderTable(NULL)
    output$leave_one_out_results <- renderTable(NULL)
    output$mr_plot <- renderPlot(NULL)
    output$funnel_plot <- renderPlot(NULL)
    output$loo_forest_plot <- renderPlot(NULL)
    showNotification("Uploads reset successfully!", type = "message")
  })
}

shinyApp(ui = ui, server = server)
