library(shiny)
library(shinyjs)
library(TwoSampleMR)
library(dplyr)
library(rmarkdown)
library(tools)

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

ui <- fluidPage(
  useShinyjs(),  # Ensure this is present so shinyjs::reset() works
  titlePanel("Two-Sample Mendelian Randomization (MR) with Sensitivity Analyses"),
  sidebarLayout(
    sidebarPanel(
      h4("Upload GWAS Summary Statistics"),
      fileInput("exposure_file", "Upload Exposure GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      fileInput("outcome_file", "Upload Outcome GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      selectInput("sep", "File Separator:",
                  choices = c("Tab" = "\t", "Comma" = ","), 
                  selected = "\t"),
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
    file_info <- input$exposure_file
    if (!is.null(file_info)) {
      ext <- tolower(file_ext(file_info$name))
      if (ext == "csv") {
        updateSelectInput(session, "sep", selected = ",")
      } else if (ext %in% c("tsv", "txt")) {
        updateSelectInput(session, "sep", selected = "\t")
      }
      rv$exposure_data <- read_uploaded_file(file_info, input$sep)
    }
  })
  
  observeEvent(input$outcome_file, {
    file_info <- input$outcome_file
    if (!is.null(file_info)) {
      ext <- tolower(file_ext(file_info$name))
      if (ext == "csv") {
        updateSelectInput(session, "sep", selected = ",")
      } else if (ext %in% c("tsv", "txt")) {
        updateSelectInput(session, "sep", selected = "\t")
      }
      rv$outcome_data <- read_uploaded_file(file_info, input$sep)
    }
  })
  
  observeEvent(input$update_threshold, {
    rv$pval_threshold <- input$pval_threshold
  })
  
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
  
  harmonized_data <- eventReactive(input$run_analysis, {
    req(rv$exposure_data, rv$outcome_data)
    exposure_filtered <- rv$exposure_data %>%
      filter(get(input$exposure_pval) < rv$pval_threshold)
    harmonise_data(exposure_dat = exposure_filtered, outcome_dat = rv$outcome_data)
  })
  
  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr(harmonized_data())
  })
  
  output$harmonized_data <- renderTable({
    req(harmonized_data())
    head(harmonized_data())
  })
  
  output$mr_results <- renderTable({
    req(mr_results())
    mr_results()
  })
  
  observeEvent(input$reset_uploads, {
    rv$exposure_data <- NULL
    rv$outcome_data <- NULL
    
    # Visually clear the file inputs
    reset("exposure_file")
    reset("outcome_file")
    
    # Clear displayed tables or outputs
    output$harmonized_data <- renderTable(NULL)
    output$mr_results <- renderTable(NULL)
    
    showNotification("Uploads reset successfully!", type = "message")
  })
}

shinyApp(ui = ui, server = server)
