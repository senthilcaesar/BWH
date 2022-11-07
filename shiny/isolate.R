library(shiny)
reactlog::reactlog_enable()

ui <- fluidPage(
  numericInput("x", "x", value = 50, min = 0, max = 100),
  actionButton("capture", "capture"),
  textOutput("out")
)

server <- function(input, output, session) {
  
  observeEvent(input$capture, ignoreInit = TRUE, {
    output$out <- renderText({paste("You have selected", isolate(input$x))})
  })

  
}

shinyApp(ui, server)
