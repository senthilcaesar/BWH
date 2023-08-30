library(shiny)
reactiveConsole(TRUE)

ui <- fluidPage(
  numericInput("a", "a", value = 10),
  numericInput("b", "b", value = 1),
  numericInput("c", "c", value = 1),
  plotOutput("x"),
  tableOutput("y"),
  textOutput("z")
)

server <- function(input, output, session) {
  
  # 1st reactive expression
  rng <- reactive(input$a * 2)
  
  # 2nd reactive expression
  # smp reactive expression read another reactive expression rng()
  smp <- reactive(sample(rng(), input$b, replace = TRUE))
  
  # 3rd reactive expression
  bc <- reactive(input$b * input$c)
  
  # Reactive outputs
  output$x <- renderPlot(hist(smp()))
  output$y <- renderTable(max(smp()))
  output$z <- renderText(bc())
}

shinyApp(ui, server)


