# https://www.stats.govt.nz/large-datasets/csv-files-for-download/

library(shiny)

ui <- fluidPage(
  textInput("url", "Enter url"),
  actionButton("do", "Download File")
)

server <- function(input, output) {
  
  output$value <- renderText({ input$url })
  
  observeEvent(input$do, {
    req(input$url)
    destfile <- paste0("/Users/sq566/tmp/",basename(input$url))
    download.file(input$url, destfile, method="libcurl")
    })
}

shinyApp(ui, server)
