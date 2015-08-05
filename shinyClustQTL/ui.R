
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)
library(markdown)
shinyUI(navbarPage("clustQTL",
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
       tags$head(HTML("<script type='text/javascript' src='http://www.biodalliance.org/release-0.11/dalliance-compiled.js'></script>")),
       tags$head(HTML("<script type='text/javascript' src='genome.js'></script>")),
       HTML("<div id='svgHolder'></div>")
       )
    )
  )
)
