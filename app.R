library(shiny)
library(shinyjs)
library(httr)
library(jsonlite)
library(TwoSampleMR)
library(dplyr)

# ChatGPT API Call Function
get_chatgpt_response <- function(user_query) {
  api_url <- "https://api.openai.com/v1/chat/completions"
  
  # Replace "your-api-key-here" with your actual API key
  api_key <- "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"  # Replace this string with your API key
  
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
  titlePanel("Mendelian Randomization with Column Mapping"),
  sidebarLayout(
    sidebarPanel(
      h4("Upload GWAS Summary Statistics Files"),
      fileInput("exposure_file", "Upload Exposure File (TSV):", accept = c(".tsv")),
      uiOutput("exposure_columns"),
      fileInput("outcome_file", "Upload Outcome File (TSV):", accept = c(".tsv")),
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
        tabPanel("Plots", plotOutput("mr_plot"))
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Read and process uploaded files
  exposure_data <- reactive({
    req(input$exposure_file)
    read.delim(input$exposure_file$datapath, header = TRUE, stringsAsFactors = FALSE)
  })
  
  outcome_data <- reactive({
    req(input$outcome_file)
    read.delim(input$outcome_file$datapath, header = TRUE, stringsAsFactors = FALSE)
  })
  
  # Generate dropdown menus for column selection
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
