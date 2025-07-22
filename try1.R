library(shiny)
library(dplyr)
library(DT)
library(stringr)
library(TwoSampleMR)
library(ggplot2)

# Define `ao` from available outcomes
ao <- available_outcomes() |> 
  tibble::as_tibble()

# Define UI
ui <- fluidPage(
  titlePanel("ShinyMR"),

  tabsetPanel(
    tabPanel("Trait Selection",
      sidebarLayout(
        sidebarPanel(
          textInput("exposure_trait", "Exposure Trait Keyword:", value = ""),
          textInput("outcome_trait", "Outcome Trait Keyword:", value = ""),
          br(),
          verbatimTextOutput("selected_ids")
        ),

        mainPanel(
          h4("Filtered Exposure Trait Rows"),
          DTOutput("exposure_table"),
          h4("Filtered Outcome Trait Rows"),
          DTOutput("outcome_table")
        )
      )
    ),
    tabPanel("Instruments",
      h4("Extracted Instruments for Exposure Trait"),
      DTOutput("instruments_table")
    ),
    tabPanel("Outcome Data",
      h4("Extracted Outcome Data for Selected SNPs"),
      DTOutput("outcome_data_table")
    ),
    tabPanel("MR Results",
      h4("Mendelian Randomization Results"),
      DTOutput("mr_results_table")
    ),
    tabPanel("Scatter Plot",
      h4("MR Scatter Plot"),
      plotOutput("scatter_plot")
    ),
    tabPanel("Download Report",
      h4("Generate MR HTML Report"),
      downloadButton("download_report", "Download MR Report")
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  filtered_exposure <- reactive({
    req(input$exposure_trait)
    ao |> 
      filter(str_detect(trait, regex(input$exposure_trait, ignore_case = TRUE)))
  })

  filtered_outcome <- reactive({
    req(input$outcome_trait)
    ao |> 
      filter(str_detect(trait, regex(input$outcome_trait, ignore_case = TRUE)))
  })

  output$exposure_table <- renderDT({
    datatable(
      filtered_exposure(), 
      selection = list(mode = "single", target = "row"),
      options = list(ordering = TRUE)
    )
  })

  output$outcome_table <- renderDT({
    datatable(
      filtered_outcome(), 
      selection = list(mode = "single", target = "row"),
      options = list(ordering = TRUE)
    )
  })

  selected_exposure_id <- reactive({
    selected <- input$exposure_table_rows_selected
    if (length(selected)) {
      filtered_exposure()[selected, ]$id
    } else {
      NA
    }
  })

  selected_outcome_id <- reactive({
    selected <- input$outcome_table_rows_selected
    if (length(selected)) {
      filtered_outcome()[selected, ]$id
    } else {
      NA
    }
  })

  output$selected_ids <- renderText({
    paste0(
      "Selected Exposure ID: ", selected_exposure_id(), "\n",
      "Selected Outcome ID: ", selected_outcome_id()
    )
  })

  instruments <- reactive({
    id <- selected_exposure_id()
    req(!is.na(id))
    extract_instruments(outcomes = id)
  })

  outcome_data <- reactive({
    snps <- instruments()$SNP
    outcome_id <- selected_outcome_id()
    req(length(snps) > 0, !is.na(outcome_id))
    extract_outcome_data(snps = snps, outcomes = outcome_id)
  })

  harmonised_data <- reactive({
    harmonise_data(exposure_dat = instruments(), outcome_dat = outcome_data())
  })

  mr_results <- reactive({
    results <- mr(harmonised_data())
    method_info <- mr_method_list()
    left_join(results, method_info, by = c("method" = "name"))
  })

  output$instruments_table <- renderDT({
    datatable(instruments())
  })

  output$outcome_data_table <- renderDT({
    datatable(outcome_data())
  })

  output$mr_results_table <- renderDT({
    datatable(mr_results())
  })

  output$scatter_plot <- renderPlot({
    plot <- mr_scatter_plot(mr_results(), harmonised_data())
    plot[[1]]
  })

  output$download_report <- downloadHandler(
    filename = function() {
      paste0("mr_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      tmp_file <- tempfile(fileext = ".html")
      mr_report(dat = harmonised_data(), output_path = tmp_file)
      file.copy(tmp_file, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
