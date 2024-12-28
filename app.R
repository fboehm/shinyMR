library(shiny)
library(shinyjs)
library(httr)
library(jsonlite)
library(TwoSampleMR)
library(dplyr)

ui <- fluidPage(
  useShinyjs(),  # Enable JavaScript
  titlePanel("Two-Sample Mendelian Randomization (MR) Analysis with Chat Support"),
  sidebarLayout(
    sidebarPanel(
      h4("MR Analysis"),
      selectInput("exposure_trait", "Select Exposure Trait from GWAS Catalog:", choices = NULL),
      selectInput("outcome_trait", "Select Outcome Trait from GWAS Catalog:", choices = NULL),
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

server <- function(input, output, session) {
  
  # Function to fetch GWAS Catalog traits
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
  
  # Update dropdown choices for traits
  observe({
    traits <- gwas_traits()
    updateSelectInput(session, "exposure_trait", choices = traits)
    updateSelectInput(session, "outcome_trait", choices = traits)
  })
  
  # Function to fetch GWAS data for a specific trait
  fetch_gwas_data <- function(trait) {
    tryCatch({
      response <- GET(paste0("https://www.ebi.ac.uk/gwas/rest/api/traits/", URLencode(trait), "/associations"))
      if (response$status_code == 200) {
        data <- fromJSON(content(response, as = "text"), flatten = TRUE)
        # Convert data to required format for MR analysis
        processed_data <- data.frame(
          snp = data$rsId,
          beta = data$beta,
          se = data$standardError,
          effect_allele = data$effectAllele,
          other_allele = data$otherAllele,
          eaf = data$eaf,
          pval = data$pValue,
          samplesize = data$sampleSize
        )
        processed_data
      } else {
        stop("Failed to fetch GWAS summary statistics.")
      }
    }, error = function(e) {
      return(NULL)
    })
  }
  
  # Reactive values for exposure and outcome data
  exposure_data <- reactive({
    req(input$exposure_trait)
    fetch_gwas_data(input$exposure_trait)
  })
  
  outcome_data <- reactive({
    req(input$outcome_trait)
    fetch_gwas_data(input$outcome_trait)
  })
  
  # Harmonize data
  harmonized_data <- eventReactive(input$run_analysis, {
    req(exposure_data(), outcome_data())
    harmonise_data(
      exposure_dat = exposure_data(),
      outcome_dat = outcome_data()
    )
  })
  
  # Run MR analysis
  mr_results <- eventReactive(input$run_analysis, {
    req(harmonized_data())
    mr(harmonized_data())
  })
  
  # Outputs for MR analysis
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
  
  # ChatGPT integration
  chat_history <- reactiveVal("")
  
  observeEvent(input$submit_query, {
    user_message <- input$user_query
    if (nchar(user_message) > 0) {
      chat_history(paste(chat_history(), "<b>You:</b> ", user_message, "<br>"))
      
      # Placeholder for ChatGPT API response (replace with actual API integration)
      assistant_reply <- paste0("You asked: ", user_message, ". This is a placeholder reply from ChatGPT.")
      
      # Update chat history
      chat_history(paste(chat_history(), "<b>Assistant:</b> ", assistant_reply, "<br>"))
      updateTextInput(session, "user_query", value = "")
    }
    
    # Update chat box
    shinyjs::html("chat_box", chat_history())
  })
}

shinyApp(ui = ui, server = server)
