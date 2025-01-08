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
  api_key <- "your-api-key-here"  # Replace this string with your API key
  
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
      fileInput("exposure_file", "Upload Exposure File", accept = c(".csv", ".tsv")),
      fileInput("outcome_file", "Upload Outcome File", accept = c(".csv", ".tsv")),
      hr(),
      h4("Map Columns (Exposure File)"),
      uiOutput("exposure_mappings"),
      hr(),
      h4("Map Columns (Outcome File)"),
      uiOutput("outcome_mappings"),
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
  # Reactive data for column names
  exposure_columns <- reactive({
    req(input$exposure_file)
    file <- input$exposure_file
    data <- read.csv(file$datapath, nrows = 1)  # Read only the first row
    colnames(data)
  })
  
  outcome_columns <- reactive({
    req(input$outcome_file)
    file <- input$outcome_file
    data <- read.csv(file$datapath, nrows = 1)  # Read only the first row
    colnames(data)
  })
  
  # Render dropdown menus for exposure file column mappings
  output$exposure_mappings <- renderUI({
    req(exposure_columns())
    tags$div(
      selectInput("exposure_snp", "SNP Column:", choices = exposure_columns(), selected = NULL),
      selectInput("exposure_beta", "Beta Column:", choices = exposure_columns(), selected = NULL),
      selectInput("exposure_se", "SE Column:", choices = exposure_columns(), selected = NULL),
      selectInput("exposure_effect_allele", "Effect Allele Column:", choices = exposure_columns(), selected = NULL),
      selectInput("exposure_other_allele", "Other Allele Column:", choices = exposure_columns(), selected = NULL),
      selectInput("exposure_eaf", "EAF Column:", choices = exposure_columns(), selected = NULL),
      selectInput("exposure_pval", "P-value Column:", choices = exposure_columns(), selected = NULL),
      selectInput("exposure_samplesize", "Sample Size Column:", choices = exposure_columns(), selected = NULL)
    )
  })
  
  # Render dropdown menus for outcome file column mappings
  output$outcome_mappings <- renderUI({
    req(outcome_columns())
    tags$div(
      selectInput("outcome_snp", "SNP Column:", choices = outcome_columns(), selected = NULL),
      selectInput("outcome_beta", "Beta Column:", choices = outcome_columns(), selected = NULL),
      selectInput("outcome_se", "SE Column:", choices = outcome_columns(), selected = NULL),
      selectInput("outcome_effect_allele", "Effect Allele Column:", choices = outcome_columns(), selected = NULL),
      selectInput("outcome_other_allele", "Other Allele Column:", choices = outcome_columns(), selected = NULL),
      selectInput("outcome_eaf", "EAF Column:", choices = outcome_columns(), selected = NULL),
      selectInput("outcome_pval", "P-value Column:", choices = outcome_columns(), selected = NULL),
      selectInput("outcome_samplesize", "Sample Size Column:", choices = outcome_columns(), selected = NULL)
    )
  })
  
  # Read and process the uploaded files based on user-selected column mappings
  process_file <- function(file, snp, beta, se, effect_allele, other_allele, eaf, pval, samplesize) {
    req(file)
    data <- read.csv(file$datapath)
    data <- data.frame(
      snp = data[[snp]],
      beta = data[[beta]],
      se = data[[se]],
      effect_allele = data[[effect_allele]],
      other_allele = data[[other_allele]],
      eaf = data[[eaf]],
      pval = data[[pval]],
      samplesize = data[[samplesize]]
    )
    return(data)
  }
  
  exposure_data <- reactive({
    req(input$exposure_file, input$exposure_snp, input$exposure_beta, input$exposure_se, 
        input$exposure_effect_allele, input$exposure_other_allele, 
        input$exposure_eaf, input$exposure_pval, input$exposure_samplesize)
    process_file(
      input$exposure_file, input$exposure_snp, input$exposure_beta, input$exposure_se,
      input$exposure_effect_allele, input$exposure_other_allele, 
      input$exposure_eaf, input$exposure_pval, input$exposure_samplesize
    )
  })
  
  outcome_data <- reactive({
    req(input$outcome_file, input$outcome_snp, input$outcome_beta, input$outcome_se, 
        input$outcome_effect_allele, input$outcome_other_allele, 
        input$outcome_eaf, input$outcome_pval, input$outcome_samplesize)
    process_file(
      input$outcome_file, input$outcome_snp, input$outcome_beta, input$outcome_se,
      input$outcome_effect_allele, input$outcome_other_allele, 
      input$outcome_eaf, input$outcome_pval, input$outcome_samplesize
    )
  })
  
  harmonized_data <- eventReactive(input$run_analysis, {
    req(exposure_data(), outcome_data())
    harmonise_data(exposure_dat = exposure_data(), outcome_dat = outcome_data())
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
