library(shiny)
library(TwoSampleMR)
library(dplyr)

ui <- fluidPage(
  titlePanel("Two-Sample Mendelian Randomization (MR) Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("exposure", "Upload Exposure GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      fileInput("outcome", "Upload Outcome GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      selectInput("sep", "File Separator:", choices = c("Tab" = "\t", "Comma" = ","), selected = "\t"),
      actionButton("run_analysis", "Run MR Analysis")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview", 
                 h3("Exposure Data"),
                 tableOutput("preview_exposure"),
                 h3("Outcome Data"),
                 tableOutput("preview_outcome")),
        tabPanel("Harmonized Data", tableOutput("harmonized_data")),
        tabPanel("MR Results", tableOutput("mr_results")),
        tabPanel("Plots", plotOutput("mr_plot")),
        tabPanel("Debugging Info", 
                 h3("Exposure Data Columns"),
                 verbatimTextOutput("exposure_cols"),
                 h3("Outcome Data Columns"),
                 verbatimTextOutput("outcome_cols"),
                 h3("Processed Exposure Data"),
                 tableOutput("debug_exposure"),
                 h3("Processed Outcome Data"),
                 tableOutput("debug_outcome"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Column mapping function
  map_columns <- function(data, column_mapping) {
    colnames(data) <- tolower(colnames(data)) # Normalize column names
    for (col in names(column_mapping)) {
      if (column_mapping[col] %in% colnames(data)) {
        colnames(data)[colnames(data) == column_mapping[col]] <- col
      }
    }
    return(data)
  }
  
  # Column mappings for exposure and outcome data
  exposure_mapping <- list(
    snp = "rsid",
    beta = "effect",
    se = "se",
    effect_allele = "a1",
    other_allele = "a2",
    eaf = "a1_freq",
    pval = "p-value",
    samplesize = "n"
  )
  
  outcome_mapping <- list(
    snp = "rsid",
    beta = "effect",
    se = "se",
    effect_allele = "a1",
    other_allele = "a2",
    eaf = "a1_freq",
    pval = "p-value",
    samplesize = "n"
  )
  
  # Debugging: Display column names of uploaded files
  output$exposure_cols <- renderPrint({
    req(input$exposure)
    tryCatch({
      data <- read.table(input$exposure$datapath, header = TRUE, sep = input$sep, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
      colnames(data)
    }, error = function(e) {
      paste("Error reading exposure data:", e$message)
    })
  })
  
  output$outcome_cols <- renderPrint({
    req(input$outcome)
    tryCatch({
      data <- read.table(input$outcome$datapath, header = TRUE, sep = input$sep, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
      colnames(data)
    }, error = function(e) {
      paste("Error reading outcome data:", e$message)
    })
  })
  
  # Handle row mismatches in Exposure Data
  processed_exposure <- reactive({
    req(input$exposure)
    tryCatch({
      data <- read.table(input$exposure$datapath, header = TRUE, sep = input$sep, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
      map_columns(data, exposure_mapping)
    }, error = function(e) {
      stop("Error processing exposure data: ", e$message)
    })
  })
  
  # Handle row mismatches in Outcome Data
  processed_outcome <- reactive({
    req(input$outcome)
    tryCatch({
      data <- read.table(input$outcome$datapath, header = TRUE, sep = input$sep, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
      map_columns(data, outcome_mapping)
    }, error = function(e) {
      stop("Error processing outcome data: ", e$message)
    })
  })
  
  # Harmonize data
  harmonized_data <- eventReactive(input$run_analysis, {
    req(processed_exposure(), processed_outcome())
    tryCatch({
      harmonised <- harmonise_data(
        exposure_dat = processed_exposure(),
        outcome_dat = processed_outcome()
      )
      if (nrow(harmonised) == 0) stop("Harmonized data contains no valid rows.")
      return(harmonised)
    }, error = function(e) {
      stop("Error during harmonization: ", e$message)
    })
  })
  
  # Run MR analysis
  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    tryCatch({
      results <- mr(harmonized_data())
      if (nrow(results) == 0) stop("MR analysis returned no results.")
      return(results)
    }, error = function(e) {
      stop("Error during MR analysis: ", e$message)
    })
  })
  
  # Outputs
  output$preview_exposure <- renderTable({
    req(processed_exposure())
    head(processed_exposure())
  })
  
  output$preview_outcome <- renderTable({
    req(processed_outcome())
    head(processed_outcome())
  })
  
  output$harmonized_data <- renderTable({
    req(harmonized_data())
    head(harmonized_data())
  })
  
  output$mr_results <- renderTable({
    req(mr_results())
    mr_results()
  })
  
  output$mr_plot <- renderPlot({
    req(mr_results(), harmonized_data())
    mr_scatter_plot(mr_results(), harmonized_data())
  })
  
  output$debug_exposure <- renderTable({
    req(processed_exposure())
    head(processed_exposure())
  })
  
  output$debug_outcome <- renderTable({
    req(processed_outcome())
    head(processed_outcome())
  })
}

shinyApp(ui = ui, server = server)
