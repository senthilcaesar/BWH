# Complete the app below with a server function that updates out 
# with the value of x only when the capture button is pressed.

library(shiny)
reactiveConsole(TRUE)

ui <- fluidPage(
  numericInput("x", "x", value = 50, min = 0, max = 100),
  actionButton("capture", "capture"),
  textOutput("out")
)

server <- function(input, output, session) {
  observeEvent(input$capture, {
    output$out <- renderText(isolate(input$x))
  })
}

shinyApp(ui, server)
