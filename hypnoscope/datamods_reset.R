library(shiny)
library(datamods)

ui <- fluidPage(
  tags$h3("Import data with copy & paste"),
  fluidRow(
    column(
      width = 4,
      import_copypaste_ui("myid")
    ),
    column(
      width = 8,
      actionButton("reset", "Reset data"),
      tags$b("Data:"),
      verbatimTextOutput(outputId = "data")
    )
  )
)

server <- function(input, output, session) {
  
  imported <- import_copypaste_server("myid", trigger_return = "change", reset = reactive(input$reset))
  output$data <- renderPrint({
    imported$data()
  })
  
}

shinyApp(ui, server)