library(shiny)
library(shinyjs)
library(httr)
library(jsonlite)
library(TwoSampleMR)
library(dplyr)

# ChatGPT API Call Function
get_chatgpt_response <- function(user_query) {
  api_url <- "https://api.openai.com/v1/chat/completions"
  
  api_key <- "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" 
  
  body <- list(
    model = "gpt-4",
    messages = list(
      list(role = "system", content = "You are a helpful assistant for a Shiny app."),
      list(role = "user", content = user_query)
    )
  )
  
  response <- POST(
    api_url,
    add_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ),
    body = toJSON(body, auto_unbox = TRUE),
    encode = "json"
  )
  
  if (response$status_code == 200) {
    content <- content(response, as = "parsed")
    return(content$choices[[1]]$message$content)  # Return ChatGPT response
  } else {
    stop(paste("API request failed (status code:", response$status_code, ")"))
  }
}

# UI Definition
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Two-Sample Mendelian Randomization Analysis with ChatGPT"),
  sidebarLayout(
    sidebarPanel(
      h4("Select GWAS Traits"),
      selectInput("exposure_trait", "Exposure Trait (GWAS Catalog):", choices = NULL),
      selectInput("outcome_trait", "Outcome Trait (GWAS Catalog):", choices = NULL),
      actionButton("fetch_data", "Fetch GWAS Data"),
      hr(),
      h4("Column Mapping for Exposure"),
      uiOutput("exposure_columns"),
      h4("Column Mapping for Outcome"),
      uiOutput("outcome_columns"),
      actionButton("run_analysis", "Run MR Analysis"),
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
        tabPanel("Plots", plotOutput("mr_plot")),
        tabPanel("Sensitivity Analyses", tableOutput("sensitivity_results"))
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Fetch GWAS Catalog traits dynamically
  gwas_traits <- reactive({
    tryCatch({
      response <- GET("https://www.ebi.ac.uk/gwas/rest/api/traits")
      if (response$status_code == 200) {
        traits <- fromJSON(content(response, as = "text"), flatten = TRUE)$trait
        unique(traits)
      } else {
        warning("GWAS Catalog API unavailable. Using fallback traits.")
        c("Trait1", "Trait2", "Trait3")  # Example fallback traits
      }
    }, error = function(e) {
      showNotification("Error fetching GWAS Catalog traits. Using fallback list.", type = "warning")
      c("Trait1", "Trait2", "Trait3")  # Example fallback traits
    })
  })
  
  # Populate GWAS trait dropdowns when the app starts
  observe({
    traits <- gwas_traits()
    updateSelectInput(session, "exposure_trait", choices = traits)
    updateSelectInput(session, "outcome_trait", choices = traits)
  })
  
  # Fetch GWAS data when the user clicks "Fetch Data"
  query_gwas_catalog <- function(trait) {
    base_url <- "https://www.ebi.ac.uk/gwas/rest/api/traits/"
    url <- paste0(base_url, URLencode(trait), "/associations")
    
    tryCatch({
      response <- GET(url)
      if (response$status_code == 200) {
        data <- fromJSON(content(response, as = "text"), flatten = TRUE)
        data.frame(
          snp = data$rsId,
          beta = data$beta,
          se = data$standardError,
          effect_allele = data$effectAllele,
          other_allele = data$otherAllele,
          eaf = data$eaf,
          pval = data$pValue,
          samplesize = data$sampleSize
        )
      } else {
        stop("Failed to fetch GWAS data from the catalog.")
      }
    }, error = function(e) {
      showNotification("Error fetching GWAS data. Please check the GWAS catalog.", type = "error")
      return(NULL)
    })
  }
  
  exposure_data <- eventReactive(input$fetch_data, {
    req(input$exposure_trait)
    query_gwas_catalog(input$exposure_trait)
  })
  
  outcome_data <- eventReactive(input$fetch_data, {
    req(input$outcome_trait)
    query_gwas_catalog(input$outcome_trait)
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
  
  # Map columns based on user selection
  map_columns <- function(data, mapping) {
    colnames(data) <- tolower(colnames(data))
    mapped_data <- data.frame(
      snp = data[[mapping$snp]],
      beta = data[[mapping$beta]],
      se = data[[mapping$se]],
      effect_allele = data[[mapping$effect_allele]],
      other_allele = data[[mapping$other_allele]],
      eaf = data[[mapping$eaf]],
      pval = data[[mapping$pval]],
      samplesize = data[[mapping$samplesize]]
    )
    mapped_data
  }
  
  # Harmonized data
  harmonized_data <- eventReactive(input$run_analysis, {
    req(exposure_data(), outcome_data())
    
    exposure_mapping <- list(
      snp = input$exposure_snp,
      beta = input$exposure_beta,
      se = input$exposure_se,
      effect_allele = input$exposure_effect_allele,
      other_allele = input$exposure_other_allele,
      eaf = input$exposure_eaf,
      pval = input$exposure_pval,
      samplesize = input$exposure_samplesize
    )
    
    outcome_mapping <- list(
      snp = input$outcome_snp,
      beta = input$outcome_beta,
      se = input$outcome_se,
      effect_allele = input$outcome_effect_allele,
      other_allele = input$outcome_other_allele,
      eaf = input$outcome_eaf,
      pval = input$outcome_pval,
      samplesize = input$outcome_samplesize
    )
    
    exposure_mapped <- map_columns(exposure_data(), exposure_mapping)
    outcome_mapped <- map_columns(outcome_data(), outcome_mapping)
    harmonise_data(exposure_dat = exposure_mapped, outcome_dat = outcome_mapped)
  })
  
  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr(harmonized_data())
  })
  
  # Outputs for MR Analysis and Plots
  output$harmonized_data <- renderTable({ req(harmonized_data()); head(harmonized_data()) })
  output$mr_results <- renderTable({ req(mr_results()); mr_results() })
  output$sensitivity_results <- renderTable({
    req(mr_results())
    # Placeholder for sensitivity analyses
    data.frame(
      Analysis = c("Leave-One-Out", "Heterogeneity Test"),
      Result = c("Placeholder Result 1", "Placeholder Result 2")
    )
  })
  output$mr_plot <- renderPlot({ req(mr_results(), harmonized_data()); mr_scatter_plot(mr_results(), harmonized_data()) })
  
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
