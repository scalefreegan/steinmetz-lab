
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)

fluidPage(
  
  title = "Poly(A) Isoforms",
  h1("Poly(A) Isoform Viewer"),
  hr(),
  fluidRow(
    column(6,DT::dataTableOutput('dt')),
    column(6,plotOutput('manhattan',
                        click = "plot_click",
                        brush = "plot_brush"))
  ),
  fluidRow(
    verbatimTextOutput("info")
    #column(12,plotOutput('pqtl'))
  ),
  fluidRow(
    column(12,plotOutput('pqtl'))
  )
)
