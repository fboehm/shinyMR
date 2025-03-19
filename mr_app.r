library(shiny)
library(shinyjs)
library(TwoSampleMR)
library(dplyr)
library(rmarkdown)
library(tools)
library(httr)
library(jsonlite)
library(ieugwasr)  # for querying OpenGWAS

# --------------------------
# Helper Functions
# --------------------------

# 1) Read uploaded file safely
read_uploaded_file <- function(file, sep) {
  req(file)
  tryCatch({
    data <- read.table(file$datapath, header = TRUE, sep = sep,
                       stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
    if (ncol(data) < 8) {
      stop("Error: The uploaded file does not have enough columns. Please check the file format.")
    }
    return(data)
  }, error = function(e) {
    showNotification("Error reading file. Please check the file format.", type = "error")
    return(NULL)
  })
}

# 2) OLS-based lookup: common trait name -> EFO ID
get_efo_from_trait_name <- function(trait_name) {
  if (trait_name == "") return(NULL)
  base_url <- "https://www.ebi.ac.uk/ols/api/search"
  query_params <- list(q = trait_name, ontology = "efo", rows = 1)
  
  out <- tryCatch({
    response <- GET(base_url, query = query_params)
    if (response$status_code != 200) {
      return(NULL)
    }
    cont <- content(response, as = "parsed", type = "application/json")
    if (!is.null(cont$response$docs) && length(cont$response$docs) > 0) {
      iri <- cont$response$docs[[1]]$iri
      efo_id <- sub("http://www.ebi.ac.uk/efo/", "", iri, fixed = TRUE)
      return(efo_id)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
  
  out
}

# 3) OLS-based lookup: EFO ID -> common trait name
get_trait_name_from_efo <- function(efo_id) {
  if (efo_id == "") return(NULL)
  
  base_url <- "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms"
  url <- paste0(base_url, "?iri=http://www.ebi.ac.uk/efo/", efo_id)
  
  out <- tryCatch({
    response <- GET(url)
    if (response$status_code != 200) {
      return(NULL)
    }
    cont <- content(response, as = "parsed", type = "application/json")
    
    if (!is.null(cont[["_embedded"]]) &&
        !is.null(cont[["_embedded"]][["terms"]]) &&
        length(cont[["_embedded"]][["terms"]]) > 0) {
      label <- cont[["_embedded"]][["terms"]][[1]]$label
      return(label)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
  
  out
}

# --- Replacing old "GWASapi" code with "ieugwasr" equivalents ---

# A) Searching for studies by trait name using gwasinfo() from ieugwasr.
search_studies_by_trait <- function(trait_name, size = 100) {
  all_studies <- gwasinfo()
  if (is.null(all_studies)) return(NULL)
  found <- all_studies[grepl(trait_name, all_studies$trait, ignore.case = TRUE), ]
  if (nrow(found) == 0) {
    return(NULL)
  }
  head(found, size)
}

# B) Getting metadata from a list of study IDs
get_study_metadata_from_ids <- function(study_ids) {
  if (length(study_ids) == 0) return(NULL)
  metadata_list <- tryCatch({
    gwasinfo(study_ids)
  }, error = function(e) {
    NULL
  })
  if (is.null(metadata_list)) return(NULL)
  metadata_list
}

# C) Using tophits() from ieugwasr to get top associations
get_exposure_snps_ieugwas <- function(gwas_id, p_upper, size) {
  res <- tryCatch({
    hits <- tophits(id = gwas_id, pval = p_upper, clump = 0)
    hits <- hits[order(hits$p), ]
    head(hits, size)
  }, error = function(e) {
    NULL
  })
  res
}

# D) For outcome, we retrieve associations for a specific set of SNPs
get_outcome_snps_ieugwas <- function(variants, gwas_id, p_upper) {
  if (length(variants) == 0) return(NULL)
  outres <- tryCatch({
    tmp <- associations(
      variants = variants, 
      id = gwas_id, 
      proxies = 0
    )
    tmp <- tmp[tmp$p <= p_upper, ]
    tmp
  }, error = function(e) {
    NULL
  })
  outres
}

# --------------------------
# UI
# --------------------------
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Two-Sample Mendelian Randomization (MR) with Sensitivity Analyses + ieugwasr"),
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
      
      # -- "ieugwasr" based SNP retrieval --
      h4("OpenGWAS: Retrieve SNPs by Study ID"),
      textInput("exposure_trait_name", "Exposure Common Trait Name (for searching):", ""),
      actionButton("lookup_exposure_efo", "Look Up Exposure EFO (from trait name)"),
      br(),
      textInput("exposure_efo", "Exposure Study ID (e.g., ieu-a-2):", ""),
      actionButton("lookup_exposure_trait", "Look Up Trait Name (from EFO)"),
      
      numericInput("exposure_api_pval", "Max P-value for Exposure SNPs:", value = 5e-8),
      numericInput("exposure_api_size", "Number of Top Exposure Associations to Retrieve:", value = 100),
      actionButton("get_exposure_snps_api", "Get Exposure SNPs from OpenGWAS"),
      hr(),
      
      textInput("outcome_trait_name", "Outcome Common Trait Name (for searching):", ""),
      actionButton("lookup_outcome_efo", "Look Up Outcome EFO (from trait name)"),
      br(),
      textInput("outcome_efo", "Outcome Study ID (e.g., ieu-a-7):", ""),
      actionButton("lookup_outcome_trait", "Look Up Trait Name (from EFO)"),
      
      numericInput("outcome_api_pval", "Max P-value for Outcome Associations:", value = 1.0),
      actionButton("get_outcome_snps_api", "Get Outcome SNPs (for the same SNP list)"),
      hr(),
      
      # NEW UI Buttons to handle multiple-studies lookup and metadata retrieval
      h4("Retrieve Multiple GWAS & Metadata via ieugwasr"),
      actionButton("search_exposure_studies", "Search Exposure Studies by Trait Name"),
      actionButton("search_outcome_studies", "Search Outcome Studies by Trait Name"),
      hr(),
      
      # -- Column Mapping for Exposure (conditional)
      conditionalPanel(
        condition = "output.showExposureMapping",
        h4("Column Mapping for Exposure"),
        uiOutput("exposure_columns")
      ),
      
      # -- Column Mapping for Outcome (conditional)
      conditionalPanel(
        condition = "output.showOutcomeMapping",
        h4("Column Mapping for Outcome"),
        uiOutput("outcome_columns")
      ),
      hr(),
      
      actionButton("run_analysis", "Run MR & Sensitivity Analyses"),
      actionButton("reset_uploads", "Reset Uploads", style = "background-color: red; color: white;"),
      downloadButton("download_report", "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Harmonized Data", tableOutput("harmonized_data")),
        tabPanel("MR Results", tableOutput("mr_results")),
        
        tabPanel("OpenGWAS Exposure SNPs", tableOutput("exposure_snps_api")),
        tabPanel("OpenGWAS Outcome SNPs", tableOutput("outcome_snps_api")),
        
        tabPanel("Filtered Exposure SNPs (QC)", tableOutput("exposure_snp_qc")),
        
        tabPanel("Heterogeneity Test", tableOutput("heterogeneity_results")),
        tabPanel("Pleiotropy Test", tableOutput("pleiotropy_results")),
        tabPanel("Single SNP Analysis", tableOutput("single_snp_results")),
        tabPanel("Leave-One-Out Analysis", tableOutput("leave_one_out_results")),
        tabPanel("Plots", 
                 plotOutput("mr_plot"), 
                 plotOutput("funnel_plot"), 
                 plotOutput("loo_forest_plot")),
        
        tabPanel("Exposure Studies Found", tableOutput("found_exposure_studies_table")),
        tabPanel("Exposure Study Metadata", tableOutput("exposure_studies_metadata_table")),
        tabPanel("Outcome Studies Found", tableOutput("found_outcome_studies_table")),
        tabPanel("Outcome Study Metadata", tableOutput("outcome_studies_metadata_table"))
      )
    )
  )
)

# --------------------------
# Server
# --------------------------
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    exposure_data = NULL,
    outcome_data = NULL,
    pval_threshold = 5e-8,
    exposure_snps_api = NULL,
    outcome_snps_api = NULL,
    found_exposure_studies = NULL,
    exposure_studies_metadata = NULL,
    found_outcome_studies = NULL,
    outcome_studies_metadata = NULL
  )
  
  # 1) Handle local file uploads
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
  
  # 2) OLS-based lookups: common trait name <-> EFO ID
  observeEvent(input$lookup_exposure_efo, {
    if (nzchar(input$exposure_trait_name)) {
      efo_id <- get_efo_from_trait_name(input$exposure_trait_name)
      if (!is.null(efo_id)) {
        updateTextInput(session, "exposure_efo", value = efo_id)
        showNotification(paste("Found EFO ID for exposure:", efo_id), type = "message")
      } else {
        showNotification("Could not find EFO for that exposure trait name.", type = "warning")
      }
    }
  })
  
  observeEvent(input$lookup_exposure_trait, {
    if (nzchar(input$exposure_efo)) {
      trait_label <- get_trait_name_from_efo(input$exposure_efo)
      if (!is.null(trait_label)) {
        updateTextInput(session, "exposure_trait_name", value = trait_label)
        showNotification(paste("Found trait name for exposure EFO:", trait_label), type = "message")
      } else {
        showNotification("Could not find trait name for that exposure EFO ID.", type = "warning")
      }
    }
  })
  
  observeEvent(input$lookup_outcome_efo, {
    if (nzchar(input$outcome_trait_name)) {
      efo_id <- get_efo_from_trait_name(input$outcome_trait_name)
      if (!is.null(efo_id)) {
        updateTextInput(session, "outcome_efo", value = efo_id)
        showNotification(paste("Found EFO ID for outcome:", efo_id), type = "message")
      } else {
        showNotification("Could not find EFO for that outcome trait name.", type = "warning")
      }
    }
  })
  
  observeEvent(input$lookup_outcome_trait, {
    if (nzchar(input$outcome_efo)) {
      trait_label <- get_trait_name_from_efo(input$outcome_efo)
      if (!is.null(trait_label)) {
        updateTextInput(session, "outcome_trait_name", value = trait_label)
        showNotification(paste("Found trait name for outcome EFO:", trait_label), type = "message")
      } else {
        showNotification("Could not find trait name for that outcome EFO ID.", type = "warning")
      }
    }
  })
  
  # 3) ieugwasr: Retrieve exposure SNPs by Study ID
  observeEvent(input$get_exposure_snps_api, {
    req(input$exposure_efo)
    api_result <- get_exposure_snps_ieugwas(
      gwas_id = input$exposure_efo,
      p_upper = input$exposure_api_pval,
      size = input$exposure_api_size
    )
    if (!is.null(api_result) && nrow(api_result) > 0) {
      rv$exposure_snps_api <- api_result
      showNotification(paste0("Fetched top exposure SNPs for study: ", input$exposure_efo), type = "message")
    } else {
      showNotification("No SNP data returned. Check your study ID or p-value threshold.", type = "warning")
    }
  })
  
  output$exposure_snps_api <- renderTable({
    req(rv$exposure_snps_api)
    head(rv$exposure_snps_api, 50)
  })
  
  # 4) ieugwasr: Retrieve outcome SNPs for the same SNP list
  observeEvent(input$get_outcome_snps_api, {
    req(rv$exposure_snps_api, input$outcome_efo)
    exp_snps <- unique(rv$exposure_snps_api$rsid)
    
    outcome_res <- get_outcome_snps_ieugwas(
      variants = exp_snps,
      gwas_id = input$outcome_efo,
      p_upper = input$outcome_api_pval
    )
    
    rv$outcome_snps_api <- outcome_res
    if (is.null(outcome_res) || nrow(outcome_res) == 0) {
      showNotification("No outcome SNP data returned for those SNPs. Check your inputs.", type = "warning")
    } else {
      showNotification(paste0("Fetched outcome SNPs for study: ", input$outcome_efo), type = "message")
    }
  })
  
  output$outcome_snps_api <- renderTable({
    req(rv$outcome_snps_api)
    head(rv$outcome_snps_api, 50)
  })
  
  # 5) Show/hide column mappings
  output$showExposureMapping <- reactive({
    !is.null(rv$exposure_data) || !is.null(rv$exposure_snps_api)
  })
  outputOptions(output, "showExposureMapping", suspendWhenHidden = FALSE)
  
  output$showOutcomeMapping <- reactive({
    !is.null(rv$outcome_data) || !is.null(rv$outcome_snps_api)
  })
  outputOptions(output, "showOutcomeMapping", suspendWhenHidden = FALSE)
  
  # 6) Column Mappings
  output$exposure_columns <- renderUI({
    colnames_source <- NULL
    if (!is.null(rv$exposure_data)) {
      colnames_source <- names(rv$exposure_data)
    } else if (!is.null(rv$exposure_snps_api)) {
      colnames_source <- names(rv$exposure_snps_api)
    }
    if (is.null(colnames_source)) {
      colnames_source <- c("N/A")
    }
    
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
    colnames_source <- NULL
    if (!is.null(rv$outcome_data)) {
      colnames_source <- names(rv$outcome_data)
    } else if (!is.null(rv$outcome_snps_api)) {
      colnames_source <- names(rv$outcome_snps_api)
    }
    if (is.null(colnames_source)) {
      colnames_source <- c("N/A")
    }
    
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
  
  # 7) Decide which data frames to use
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
  
  # 8) Harmonize and run MR after user clicks the button
  exposure_filtered_data <- eventReactive(input$run_analysis, {
    exp_data <- final_exposure_data()
    req(exp_data)
    exp_pval_col <- input$exposure_pval
    exposure_filtered <- exp_data %>% filter(.data[[exp_pval_col]] < rv$pval_threshold)
    exposure_filtered
  })
  
  output$exposure_snp_qc <- renderTable({
    req(exposure_filtered_data())
    head(exposure_filtered_data(), 50)
  })
  
  harmonized_data <- eventReactive(input$run_analysis, {
    exp_filtered <- exposure_filtered_data()
    out_data <- final_outcome_data()
    req(exp_filtered, out_data)
    harmonise_data(
      exposure_dat = exp_filtered,
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
  
  # 9) Reset file uploads and data
  observeEvent(input$reset_uploads, {
    rv$exposure_data <- NULL
    rv$outcome_data <- NULL
    rv$exposure_snps_api <- NULL
    rv$outcome_snps_api <- NULL
    rv$found_exposure_studies <- NULL
    rv$exposure_studies_metadata <- NULL
    rv$found_outcome_studies <- NULL
    rv$outcome_studies_metadata <- NULL
    
    reset("exposure_file")
    reset("outcome_file")
    
    output$harmonized_data <- renderTable(NULL)
    output$mr_results <- renderTable(NULL)
    output$exposure_snps_api <- renderTable(NULL)
    output$outcome_snps_api <- renderTable(NULL)
    output$exposure_snp_qc <- renderTable(NULL)
    
    output$found_exposure_studies_table <- renderTable(NULL)
    output$exposure_studies_metadata_table <- renderTable(NULL)
    output$found_outcome_studies_table <- renderTable(NULL)
    output$outcome_studies_metadata_table <- renderTable(NULL)
    
    showNotification("Uploads reset successfully!", type = "message")
  })
  
  # 10) Download report
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("MR_analysis_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      # CHANGED: use "report_template.Rmd" instead of "report.Rmd"
      tempReport <- file.path(tempdir(), "report_template.Rmd")
      file.copy("report_template.Rmd", tempReport, overwrite = TRUE)
      
      params <- list(
        exposure_efo = input$exposure_efo,
        outcome_efo = input$outcome_efo,
        exposure_trait_name = input$exposure_trait_name,
        outcome_trait_name = input$outcome_trait_name,
        exposure_snps_api = rv$exposure_snps_api,
        outcome_snps_api = rv$outcome_snps_api,
        exposure_filtered_data = exposure_filtered_data(),
        harmonized_data = harmonized_data(),
        mr_results = mr_results()
      )
      
      rmarkdown::render(
        tempReport,
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv())
      )
    }
  )
  
  # NEW: Search for multiple GWAS/studies for the input trait name (Exposure)
  observeEvent(input$search_exposure_studies, {
    req(input$exposure_trait_name)
    trait_name <- input$exposure_trait_name
    
    studies_found <- search_studies_by_trait(trait_name, size = 100)
    if (is.null(studies_found) || nrow(studies_found) == 0) {
      showNotification("No exposure studies found for the given trait name in OpenGWAS.", type = "warning")
      return(NULL)
    }
    
    rv$found_exposure_studies <- studies_found
    all_ids <- studies_found$id
    meta_df <- get_study_metadata_from_ids(all_ids)
    
    if (is.null(meta_df)) {
      showNotification("No metadata could be retrieved for these exposure studies.", type = "warning")
    }
    rv$exposure_studies_metadata <- meta_df
    
    showNotification("Exposure studies and metadata retrieved successfully (via ieugwasr).", type = "message")
  })
  
  # NEW: Search for multiple GWAS/studies for the input trait name (Outcome)
  observeEvent(input$search_outcome_studies, {
    req(input$outcome_trait_name)
    trait_name <- input$outcome_trait_name
    
    studies_found <- search_studies_by_trait(trait_name, size = 100)
    if (is.null(studies_found) || nrow(studies_found) == 0) {
      showNotification("No outcome studies found for the given trait name in OpenGWAS.", type = "warning")
      return(NULL)
    }
    
    rv$found_outcome_studies <- studies_found
    all_ids <- studies_found$id
    meta_df <- get_study_metadata_from_ids(all_ids)
    
    if (is.null(meta_df)) {
      showNotification("No metadata could be retrieved for these outcome studies.", type = "warning")
    }
    rv$outcome_studies_metadata <- meta_df
    
    showNotification("Outcome studies and metadata retrieved successfully (via ieugwasr).", type = "message")
  })
  
  output$found_exposure_studies_table <- renderTable({
    req(rv$found_exposure_studies)
    rv$found_exposure_studies
  })
  
  output$exposure_studies_metadata_table <- renderTable({
    req(rv$exposure_studies_metadata)
    rv$exposure_studies_metadata
  })
  
  output$found_outcome_studies_table <- renderTable({
    req(rv$found_outcome_studies)
    rv$found_outcome_studies
  })
  
  output$outcome_studies_metadata_table <- renderTable({
    req(rv$outcome_studies_metadata)
    rv$outcome_studies_metadata
  })
}

shinyApp(ui = ui, server = server)
