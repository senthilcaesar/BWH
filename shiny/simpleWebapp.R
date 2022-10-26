sales <- readr::read_csv('/Users/sq566/myapp/neiss/sales_data_sample.csv',
                         show_col_types = FALSE, na="")
sales <- sales[c(
  "TERRITORY", "ORDERDATE", "ORDERNUMBER", "PRODUCTCODE",
  "QUANTITYORDERED", "PRICEEACH"
)]

library(shiny)

ui <- fluidPage(
  selectInput("territory", "territory", choices = unique(sales$TERRITORY)),
  tableOutput("selected")
)
server <- function(input, output, session) {
  selected <- reactive(subset(sales, TERRITORY == input$territory))
  output$selected <- renderTable(head(selected(), 10))
}

shinyApp(ui, server)