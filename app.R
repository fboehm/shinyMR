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
                 h3("Original Exposure Data Columns"),
                 verbatimTextOutput("exposure_cols"),
                 h3("Original Outcome Data Columns"),
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
  
  # Function to read data dynamically
  read_data <- function(file, sep) {
    if (is.null(file)) return(NULL)
    ext <- tools::file_ext(file$name)
    if (ext %in% c("tsv", "csv")) {
      tryCatch(
        read.table(file$datapath, header = TRUE, sep = sep, stringsAsFactors = FALSE, fill = TRUE),
        error = function(e) stop("Error reading file: ", e$message)
      )
    } else {
      stop("Unsupported file format. Please upload a TSV or CSV file.")
    }
  }
  
  # Function to standardize column names and fill missing ones
  standardize_data <- function(data, role) {
    required_cols <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval")
    colnames(data) <- tolower(gsub("\\.", "_", colnames(data)))
    
    # Ensure all required columns are present
    for (col in required_cols) {
      if (!(col %in% colnames(data))) {
        data[[col]] <- NA
      }
    }
    
    standardized_data <- data %>% select(all_of(required_cols))
    colnames(standardized_data) <- paste0(required_cols, ".", role)
    
    # Add missing mandatory columns for harmonise_data
    mandatory_cols <- c("SNP", "id", "exposure")
    for (col in mandatory_cols) {
      full_col_name <- paste0(col, ".", role)
      if (!(full_col_name %in% colnames(standardized_data))) {
        if (col == "SNP") {
          standardized_data[[full_col_name]] <- paste0("SNP_", seq_len(nrow(standardized_data)))
        } else {
          standardized_data[[full_col_name]] <- role
        }
      }
    }
    
    return(standardized_data)
  }
  
  # Reactive expressions for reading files
  exposure_data <- reactive({
    req(input$exposure)
    read_data(input$exposure, input$sep)
  })
  
  outcome_data <- reactive({
    req(input$outcome)
    read_data(input$outcome, input$sep)
  })
  
  # Harmonizing data
  harmonized_data <- eventReactive(input$run_analysis, {
    req(exposure_data(), outcome_data())
    exp_data <- standardize_data(exposure_data(), "exposure")
    out_data <- standardize_data(outcome_data(), "outcome")
    
    tryCatch({
      harmonized <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)
      if (nrow(harmonized) == 0) {
        stop("Harmonized data contains no valid rows. Check if your datasets align correctly.")
      }
      return(harmonized)
    }, error = function(e) {
      stop("Error during data harmonization: ", e$message)
    })
  })
  
  # Running MR analysis
  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    tryCatch({
      results <- mr(harmonized_data())
      if (nrow(results) == 0) {
        stop("MR analysis returned no results. Check your harmonized data.")
      }
      return(results)
    }, error = function(e) {
      stop("Error during MR analysis: ", e$message)
    })
  })
  
  # Outputs
  output$preview_exposure <- renderTable({
    req(exposure_data())
    head(exposure_data())
  })
  
  output$preview_outcome <- renderTable({
    req(outcome_data())
    head(outcome_data())
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
    req(exposure_data())
    head(standardize_data(exposure_data(), "exposure"))
  })
  
  output$debug_outcome <- renderTable({
    req(outcome_data())
    head(standardize_data(outcome_data(), "outcome"))
  })
  
  output$exposure_cols <- renderPrint({
    req(exposure_data())
    colnames(exposure_data())
  })
  
  output$outcome_cols <- renderPrint({
    req(outcome_data())
    colnames(outcome_data())
  })
}

shinyApp(ui = ui, server = server)
