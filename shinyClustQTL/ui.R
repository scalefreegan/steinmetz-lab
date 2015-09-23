
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)
library(markdown)
shinyUI(navbarPage("clustQTL", selected = "Inspect",            
         tabPanel("Detect", id = "detect", fluidPage()),
         tabPanel("Inspect",id = "inspect", fluidPage(
            title = "Poly(A) Isoforms",
            #h1("Poly(A) Isoform Viewer"),
            hr(),
            fluidRow(
              column(6,
                fluidRow(
                  column(12,sliderInput("alpha", " QTL Threshold: -log10(pval),", max = 30, min = 0, value = 5), offset = 3)
                ),
                fluidRow(
                  column(12,DT::dataTableOutput('dt'))
                )
              ),
              column(6,
                 fluidRow(
                   column(12,h4("Manhattan Plot"),offset = 4)
                 ),
                 fluidRow(
                   column(12, plotOutput('manhattan',
                              click = "plot_click",
                              brush = "plot_brush",
                              dblclick = "plot_dblclick"))
               )
            ),
            fluidRow(
              column(4,h4("Phenotype at selected location"),offset = 4)
            ),
            fluidRow(
              column(12,plotOutput('pqtl'))
              )
            )
          )
         ),
       tabPanel("Genome Browser", id = "genome", fluidPage(
         tags$head(HTML("<title>Yeast JBrowse </title>")),
         htmlOutput("link", inline = TRUE)
         #htmlOutput('link')
         #HTML("<div style='width: 100%; height: 600px'><iframe style='border: 1px solid black' src='http://localhost/~brooks/JBrowse-1.11.6/' width='100%' height='100%'></iframe></div>")
         )
      )
    )
)
