# --------------------------
# Updated Shiny App with GWASapi Integration
# --------------------------
library(shiny)
library(shinyjs)
library(TwoSampleMR)
library(dplyr)
library(rmarkdown)
library(tools)
library(GWASapi)  # New: For querying GWAS Catalog

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
  useShinyjs(),
  titlePanel("Two-Sample Mendelian Randomization (MR) with Sensitivity Analyses + GWASapi"),
  sidebarLayout(
    sidebarPanel(
      # -- Existing file upload section --
      h4("Upload GWAS Summary Statistics (Local Files)"),
      fileInput("exposure_file", "Upload Exposure GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      fileInput("outcome_file", "Upload Outcome GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      selectInput("sep", "File Separator:",
                  choices = c("Tab" = "\t", "Comma" = ","), 
                  selected = "\t"),
      numericInput("pval_threshold", "P-value threshold for SNP selection:", value = 5e-8, step = 1e-8),
      actionButton("update_threshold", "Update P-Value Threshold"),
      hr(),
      
      # -- New: GWASapi-based SNP retrieval section --
      h4("GWASapi: Retrieve SNPs by Trait"),
      textInput("exposure_efo", "Exposure Trait EFO ID (e.g., EFO_0001360):", ""),
      numericInput("exposure_api_pval", "Max P-value for Exposure SNPs:", value = 5e-8),
      numericInput("exposure_api_size", "Number of Top Exposure Associations to Retrieve:", value = 100),
      actionButton("get_exposure_snps_api", "Get Exposure SNPs from GWAS Catalog"),
      br(),
      
      textInput("outcome_efo", "Outcome Trait EFO ID (e.g., EFO_0001360):", ""),
      numericInput("outcome_api_pval", "Max P-value for Outcome Associations:", value = 1.0),
      actionButton("get_outcome_snps_api", "Get Outcome SNPs (for the same SNP list)"),
      hr(),
      
      # -- Column mapping from local or API-based data --
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
        
        # -- New: Show results from GWASapi queries --
        tabPanel("GWAS Catalog Exposure SNPs", tableOutput("exposure_snps_api")),
        tabPanel("GWAS Catalog Outcome SNPs", tableOutput("outcome_snps_api")),
        
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
  rv <- reactiveValues(
    exposure_data = NULL,      # For local file or processed data
    outcome_data = NULL,       # For local file or processed data
    pval_threshold = 5e-8,
    # New: For GWASapi queries
    exposure_snps_api = NULL,  # Exposure trait associations from GWAS Catalog
    outcome_snps_api = NULL    # Outcome trait associations for the same SNPs
  )
  
  # 1. Handle local file uploads -----------------------------------------
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
  
  
  # 2. GWASapi: Retrieve exposure SNPs by trait ---------------------------
  observeEvent(input$get_exposure_snps_api, {
    req(input$exposure_efo)
    
    # Query the GWAS Catalog for associations with a given trait EFO
    # restricting to p-value < input$exposure_api_pval
    # and returning up to input$exposure_api_size associations.
    
    # get_trait_asso() can retrieve all associations for a given trait EFO code
    # We'll filter them by p_value < input$exposure_api_pval
    # Then keep only the top 'exposure_api_size' associations by significance.
    
    api_result <- get_trait_asso(
      input$exposure_efo,
      p_upper = input$exposure_api_pval,
      size = input$exposure_api_size
    )
    
    if (!is.null(api_result) && nrow(api_result) > 0) {
      # Sort by p-value ascending
      api_result <- api_result[order(api_result$p_value), ]
    }
    
    rv$exposure_snps_api <- api_result
  })
  
  # Display the resulting exposure SNPs from the GWAS Catalog
  output$exposure_snps_api <- renderTable({
    req(rv$exposure_snps_api)
    head(rv$exposure_snps_api, 50)
  })
  
  
  # 3. GWASapi: Retrieve outcome SNPs for the same SNP list ----------------
  observeEvent(input$get_outcome_snps_api, {
    req(rv$exposure_snps_api, input$outcome_efo)
    
    # The approach here is:
    #  - Extract the distinct variant identifiers (rsID) from the exposure results
    #  - For each variant, call get_variant() from GWASapi
    #  - Filter to keep only associations matching the user-specified outcome EFO
    #  - Restrict to p_value < input$outcome_api_pval
    
    exp_snps <- unique(rv$exposure_snps_api$variant_id)  # Typically "rsID"
    results_list <- list()
    
    # A small helper to safely query variants
    query_variant_for_outcome <- function(rs) {
      # We can query all associations for the variant, then filter by outcome trait
      vres <- get_variant(rs, size = 1000)  # can adjust size if needed
      if (!is.null(vres) && nrow(vres) > 0) {
        # Filter by trait and p-value
        vres <- vres[vres$trait == input$outcome_efo & vres$p_value <= input$outcome_api_pval, ]
      }
      vres
    }
    
    # Loop over the SNPs and combine
    for (snp in exp_snps) {
      vdat <- query_variant_for_outcome(snp)
      if (!is.null(vdat) && nrow(vdat) > 0) {
        results_list[[snp]] <- vdat
      }
    }
    
    outcome_res <- if (length(results_list) > 0) {
      do.call(rbind, results_list)
    } else {
      NULL
    }
    rv$outcome_snps_api <- outcome_res
  })
  
  # Display the resulting outcome associations
  output$outcome_snps_api <- renderTable({
    req(rv$outcome_snps_api)
    head(rv$outcome_snps_api, 50)
  })
  
  
  # 4. UI outputs for column mapping ---------------------------------------
  output$exposure_columns <- renderUI({
    # We choose which data to map from: local file or possibly from the
    # rv$exposure_snps_api. We'll default to local file if that was provided,
    # otherwise let the user map from the GWASapi results, if desired.
    
    # If the user wants to do MR from the local file, we keep the old approach.
    # They could also do a separate approach: convert the API results to the required format.
    
    colnames_source <- NULL
    
    if (!is.null(rv$exposure_data)) {
      colnames_source <- names(rv$exposure_data)
    } else if (!is.null(rv$exposure_snps_api)) {
      colnames_source <- names(rv$exposure_snps_api)
    }
    
    req(colnames_source)
    
    tagList(
      selectInput("exposure_snp", "SNP Column:", choices = colnames_source),
      selectInput("exposure_beta", "Beta Column:", choices = colnames_source),
      selectInput("exposure_se", "SE Column:", choices = colnames_source),
      selectInput("exposure_effect_allele", "Effect Allele Column:", choices = colnames_source),
      selectInput("exposure_other_allele", "Other Allele Column:", choices = colnames_source),
      selectInput("exposure_eaf", "EAF Column:", choices = colnames_source),
      selectInput("exposure_pval", "P-Value Column:", choices = colnames_source),
      selectInput("exposure_samplesize", "Sample Size Column:", choices = colnames_source)
    )
  })
  
  output$outcome_columns <- renderUI({
    # Similar logic as exposure columns
    colnames_source <- NULL
    
    if (!is.null(rv$outcome_data)) {
      colnames_source <- names(rv$outcome_data)
    } else if (!is.null(rv$outcome_snps_api)) {
      colnames_source <- names(rv$outcome_snps_api)
    }
    
    req(colnames_source)
    
    tagList(
      selectInput("outcome_snp", "SNP Column:", choices = colnames_source),
      selectInput("outcome_beta", "Beta Column:", choices = colnames_source),
      selectInput("outcome_se", "SE Column:", choices = colnames_source),
      selectInput("outcome_effect_allele", "Effect Allele Column:", choices = colnames_source),
      selectInput("outcome_other_allele", "Other Allele Column:", choices = colnames_source),
      selectInput("outcome_eaf", "EAF Column:", choices = colnames_source),
      selectInput("outcome_pval", "P-Value Column:", choices = colnames_source),
      selectInput("outcome_samplesize", "Sample Size Column:", choices = colnames_source)
    )
  })
  
  
  # 5. MR pipeline --------------------------------------------------------
  # Decide which data frames to use for exposure/outcome: local or from API
  # For simplicity, we default to local if available, otherwise from the API.
  
  final_exposure_data <- reactive({
    if (!is.null(rv$exposure_data)) {
      rv$exposure_data
    } else {
      rv$exposure_snps_api
    }
  })
  
  final_outcome_data <- reactive({
    if (!is.null(rv$outcome_data)) {
      rv$outcome_data
    } else {
      rv$outcome_snps_api
    }
  })
  
  # Harmonize and run MR after user clicks "Run MR & Sensitivity Analyses"
  harmonized_data <- eventReactive(input$run_analysis, {
    exp_data <- final_exposure_data()
    out_data <- final_outcome_data()
    req(exp_data, out_data)
    
    # Filter exposure data by user-chosen p-value threshold
    # (local or from API). The input$exposure_pval is the name of the column,
    # while rv$pval_threshold is the numeric threshold.
    exposure_filtered <- exp_data %>%
      filter(get(input$exposure_pval) < rv$pval_threshold)
    
    harmonise_data(
      exposure_dat = exposure_filtered,
      outcome_dat  = out_data
    )
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
  
  
  # 6. Reset file uploads and data ----------------------------------------
  observeEvent(input$reset_uploads, {
    rv$exposure_data <- NULL
    rv$outcome_data <- NULL
    rv$exposure_snps_api <- NULL
    rv$outcome_snps_api <- NULL
    
    # Visually clear the file inputs
    reset("exposure_file")
    reset("outcome_file")
    
    # Clear displayed tables or outputs
    output$harmonized_data <- renderTable(NULL)
    output$mr_results <- renderTable(NULL)
    output$exposure_snps_api <- renderTable(NULL)
    output$outcome_snps_api <- renderTable(NULL)
    
    showNotification("Uploads reset successfully!", type = "message")
  })
}

shinyApp(ui = ui, server = server)
