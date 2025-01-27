library(shiny)
library(shinyjs)
library(httr)
library(jsonlite)
library(TwoSampleMR)
library(dplyr)
library(rmarkdown)

# Helper function: Fetch EFO traits
fetch_efo_traits <- function() {
  base_url <- "https://www.ebi.ac.uk/gwas/rest/api/traits"
  response <- GET(base_url)
  
  if (response$status_code == 200) {
    traits <- fromJSON(content(response, as = "text"), flatten = TRUE)
    return(unique(traits$trait))
  } else {
    return(NULL)
  }
}

# Helper function: Fetch associations for a given trait
fetch_associations <- function(trait) {
  base_url <- "https://www.ebi.ac.uk/gwas/rest/api/traits/"
  url <- paste0(base_url, URLencode(trait), "/associations")
  response <- GET(url)
  
  if (response$status_code == 200) {
    data <- fromJSON(content(response, as = "text"), flatten = TRUE)
    if (!is.null(data$rsId)) {
      return(data.frame(
        snp = data$rsId,
        beta = data$beta,
        se = data$standardError,
        effect_allele = data$effectAllele,
        other_allele = data$otherAllele,
        eaf = data$eaf,
        pval = data$pValue,
        samplesize = data$sampleSize
      ))
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

# UI Definition
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Searchable GWAS Traits and MR Analysis"),
  sidebarLayout(
    sidebarPanel(
      h4("Search and Select GWAS Traits"),
      selectizeInput(
        "exposure_trait",
        "Exposure Trait:",
        choices = NULL,
        options = list(
          placeholder = "Type to search for a trait",
          maxOptions = 10  # Limits the number of autocomplete options shown
        )
      ),
      selectizeInput(
        "outcome_trait",
        "Outcome Trait:",
        choices = NULL,
        options = list(
          placeholder = "Type to search for a trait",
          maxOptions = 10
        )
      ),
      actionButton("fetch_data", "Fetch GWAS Data"),
      hr(),
      h4("Column Mapping for Exposure"),
      uiOutput("exposure_columns"),
      h4("Column Mapping for Outcome"),
      uiOutput("outcome_columns"),
      actionButton("run_analysis", "Run MR Analysis"),
      downloadButton("download_report", "Download RMarkdown Report"),
      hr(),
      h4("Chat with Assistant"),
      textInput("user_query", "Ask a question:", ""),
      actionButton("submit_query", "Send"),
      div(id = "chat_box", style = "height: 200px; overflow-y: auto; border: 1px solid #ddd; padding: 10px;")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Harmonized Data", tableOutput("harmonized_data")),
        tabPanel("MR Results", tableOutput("mr_results")),
        tabPanel("Plots", plotOutput("mr_plot"))
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Dynamically fetch and populate traits
  observe({
    traits <- fetch_efo_traits()
    if (!is.null(traits)) {
      updateSelectizeInput(session, "exposure_trait", choices = traits, server = TRUE)
      updateSelectizeInput(session, "outcome_trait", choices = traits, server = TRUE)
    } else {
      showNotification("Failed to fetch GWAS traits. Using fallback list.", type = "warning")
      updateSelectizeInput(session, "exposure_trait", choices = c("Trait1", "Trait2", "Trait3"), server = TRUE)
      updateSelectizeInput(session, "outcome_trait", choices = c("Trait1", "Trait2", "Trait3"), server = TRUE)
    }
  })
  
  # Fetch GWAS data for selected traits
  exposure_data <- eventReactive(input$fetch_data, {
    req(input$exposure_trait)
    fetch_associations(input$exposure_trait)
  })
  
  outcome_data <- eventReactive(input$fetch_data, {
    req(input$outcome_trait)
    fetch_associations(input$outcome_trait)
  })
  
  # Generate dropdown menus for column mapping
  output$exposure_columns <- renderUI({
    req(exposure_data())
    column_names <- names(exposure_data())
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
    req(outcome_data())
    column_names <- names(outcome_data())
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
  
  # Harmonize and analyze data
  harmonized_data <- eventReactive(input$run_analysis, {
    req(exposure_data(), outcome_data())
    harmonise_data(exposure_dat = exposure_data(), outcome_dat = outcome_data())
  })
  
  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr(harmonized_data())
  })
  
  output$harmonized_data <- renderTable({ req(harmonized_data()); head(harmonized_data()) })
  output$mr_results <- renderTable({ req(mr_results()); mr_results() })
  output$mr_plot <- renderPlot({ req(mr_results(), harmonized_data()); mr_scatter_plot(mr_results(), harmonized_data()) })
  
  # Generate RMarkdown Report
  output$download_report <- downloadHandler(
    filename = function() { paste("MR_Analysis_Report_", Sys.Date(), ".html", sep = "") },
    content = function(file) {
      params <- list(
        harmonized_data = harmonized_data(),
        mr_results = mr_results()
      )
      rmarkdown::render("report_template.Rmd", output_file = file, params = params, envir = new.env(parent = globalenv()))
    }
  )
  
  # ChatGPT Integration
  chat_history <- reactiveVal("")
  observeEvent(input$submit_query, {
    user_message <- input$user_query
    if (nchar(user_message) > 0) {
      chat_history(paste(chat_history(), "<b>You:</b> ", user_message, "<br>"))
      tryCatch({
        assistant_reply <- get_chatgpt_response(user_message)
        chat_history(paste(chat_history(), "<b>Assistant:</b> ", assistant_reply, "<br>"))
      }, error = function(e) {
        chat_history(paste(chat_history(), "<b>Assistant:</b> ", "Error: Unable to fetch response.", "<br>"))
      })
      updateTextInput(session, "user_query", value = "")
    }
    shinyjs::html("chat_box", chat_history())
  })
}

shinyApp(ui = ui, server = server)
