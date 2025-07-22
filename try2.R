library(shiny)
library(dplyr)
library(DT)
library(stringr)
library(TwoSampleMR)
library(glue)

ui <- fluidPage(
  titlePanel("Two-sample Mendelian Randomization (MR)"),

  sidebarLayout(
    sidebarPanel(
      textInput("exp_keyword", "Exposure Trait Keyword:"),
      textInput("out_keyword", "Outcome Trait Keyword:"),
      downloadButton("download_qmd", "Download Analysis Code")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Traits",
          DTOutput("exp_table"),
          DTOutput("out_table")
        ),
        tabPanel("Exposure Data", DTOutput("exp_data_table")),
        tabPanel("Outcome Data", DTOutput("out_data_table")),
        tabPanel("MR Results", DTOutput("mr_results")),
        tabPanel("Plot", plotOutput("mr_plot")),
        tabPanel("Sensitivity Analyses",
                 plotOutput("funnel_plot"),
                 plotOutput("leaveoneout_plot"))
      )
    )
  )
)

server <- function(input, output, session) {

  ao <- available_outcomes() |> tibble::as_tibble()

  exp_data <- reactive({
    ao |> filter(str_detect(trait, regex(input$exp_keyword, ignore_case = TRUE)))
  })

  out_data <- reactive({
    ao |> filter(str_detect(trait, regex(input$out_keyword, ignore_case = TRUE)))
  })

  output$exp_table <- renderDT({ datatable(exp_data(), selection = "single") })
  output$out_table <- renderDT({ datatable(out_data(), selection = "single") })

  selected_exp <- reactive({ req(input$exp_table_rows_selected); exp_data()[input$exp_table_rows_selected,] })
  selected_out <- reactive({ req(input$out_table_rows_selected); out_data()[input$out_table_rows_selected,] })

  instruments <- reactive({ extract_instruments(outcomes = selected_exp()$id) })
  outcome_dat <- reactive({ extract_outcome_data(snps = instruments()$SNP, outcomes = selected_out()$id) })
  harmonised <- reactive({ harmonise_data(instruments(), outcome_dat()) })

  mr_res <- reactive({ mr(harmonised()) })

  output$exp_data_table <- renderDT({ datatable(instruments()) })
  output$out_data_table <- renderDT({ datatable(outcome_dat()) })

  output$mr_results <- renderDT({ datatable(mr_res()) })

  output$mr_plot <- renderPlot({ mr_scatter_plot(mr_res(), harmonised())[[1]] })
  output$funnel_plot <- renderPlot({ mr_funnel_plot(mr_singlesnp(harmonised()))[[1]] })
  output$leaveoneout_plot <- renderPlot({ mr_leaveoneout_plot(mr_leaveoneout(harmonised()))[[1]] })

  user_code <- reactive({
    glue::glue(
      """
      ---
      title: "Two-sample MR Analysis"
      format: html
      ---

      ```{{r}}
      library(TwoSampleMR)

      # Selected traits
      exp_id <- "{selected_exp()$id}"
      out_id <- "{selected_out()$id}"

      # Extract instruments and outcome data
      exposure <- extract_instruments(outcomes = exp_id)
      outcome <- extract_outcome_data(snps = exposure$SNP, outcomes = out_id)

      # Harmonize datasets
      harmonised <- harmonise_data(exposure, outcome)

      # Perform MR analysis
      mr_results <- mr(harmonised)
      mr_results

      # Plot MR results
      mr_scatter_plot(mr_results, harmonised)

      # Sensitivity analyses
      singlesnp_results <- mr_singlesnp(harmonised)
      mr_funnel_plot(singlesnp_results)

      loo_results <- mr_leaveoneout(harmonised)
      mr_leaveoneout_plot(loo_results)
      ```
      """
    )
  })

  output$download_qmd <- downloadHandler(
    filename = function() paste0("mr_analysis_", Sys.Date(), ".qmd"),
    content = function(file) {
      writeLines(user_code(), con = file)
    }
  )
}

shinyApp(ui, server)
