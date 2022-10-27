library(ggplot2)
library(shiny)

ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  textInput("dataset", "Dataset name"), 
  tableOutput("data")
)

server <- function(input, output, session) {
  data <- reactive({
    req(input$dataset)
    
    exists <- exists(input$dataset, "package:datasets")
    shinyFeedback::feedbackDanger("dataset", !exists, "Unknown dataset")
    req(exists, cancelOutput = TRUE)
    
    get(input$dataset, "package:datasets")
  })
  
  output$data <- renderTable({
    head(data())
  })
}

shinyApp(ui, server)
