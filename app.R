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
  titlePanel("Mendelian Randomization with Chat Support"),
  sidebarLayout(
    sidebarPanel(
      h4("MR Analysis"),
      selectInput("exposure_trait", "Select Exposure Trait:", choices = NULL),
      selectInput("outcome_trait", "Select Outcome Trait:", choices = NULL),
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
  # Fetch GWAS Catalog Traits Dynamically
  gwas_traits <- reactive({
    tryCatch({
      response <- GET("https://www.ebi.ac.uk/gwas/rest/api/associations")
      if (response$status_code == 200) {
        traits <- fromJSON(content(response, as = "text"), flatten = TRUE)$trait
        unique(traits)
      } else {
        stop("Failed to fetch GWAS Catalog traits.")
      }
    }, error = function(e) {
      return(c("Error fetching traits" = NA))
    })
  })
  
  # Update Trait Selection Dropdowns
  observe({
    traits <- gwas_traits()
    updateSelectInput(session, "exposure_trait", choices = traits)
    updateSelectInput(session, "outcome_trait", choices = traits)
  })
  
  # Fetch GWAS Data for a Specific Trait
  fetch_gwas_data <- function(trait) {
    tryCatch({
      response <- GET(paste0("https://www.ebi.ac.uk/gwas/rest/api/traits/", URLencode(trait), "/associations"))
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
        stop("Failed to fetch GWAS data.")
      }
    }, error = function(e) {
      return(NULL)
    })
  }
  
  # Reactive Data Processing
  exposure_data <- reactive({ req(input$exposure_trait); fetch_gwas_data(input$exposure_trait) })
  outcome_data <- reactive({ req(input$outcome_trait); fetch_gwas_data(input$outcome_trait) })
  
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
