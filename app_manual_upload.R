library(shiny)
library(shinyjs)
library(TwoSampleMR)
library(dplyr)
library(rmarkdown)

# UI Definition
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Manual Upload of GWAS Summary Statistics for MR Analysis"),
  sidebarLayout(
    sidebarPanel(
      h4("Upload GWAS Summary Statistics"),
      fileInput("exposure_file", "Upload Exposure GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      fileInput("outcome_file", "Upload Outcome GWAS Summary Stats (TSV or CSV):", accept = c(".tsv", ".csv")),
      selectInput("sep", "File Separator:", choices = c("Tab" = "\t", "Comma" = ","), selected = "\t"),
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
  # Read uploaded files
  read_uploaded_file <- function(file, sep) {
    req(file)
    read.table(file$datapath, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  }
  
  exposure_data <- reactive({
    req(input$exposure_file)
    read_uploaded_file(input$exposure_file, input$sep)
  })
  
  outcome_data <- reactive({
    req(input$outcome_file)
    read_uploaded_file(input$outcome_file, input$sep)
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
