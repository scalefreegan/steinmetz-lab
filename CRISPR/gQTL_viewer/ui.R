
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)
library(markdown)
shinyUI(navbarPage("gQTL Explorer", selected = "Table Viewer",
         tabPanel("Table Viewer",id = "inspect", fluidPage(
            title = "gQTL Explorer",
            #h1("Poly(A) Isoform Viewer"),
            hr(),
            fluidRow(
              column(12,
                 fluidRow(
                   column(12, h2("gQTLs From Gagneur et al 2013",align = "center"))
                   ),
                 fluidRow(
                   column(12, htmlOutput("image"))
               )
              )
            ),
            fluidRow(
              column(12,
                     fluidRow(
                       column(12,sliderInput("cr", "Custom Interval (bp)", max = 100000, min = 0, value = 5, step = 100, width = "100%"))
                     )
              )
            ),
            fluidRow(
              column(1),
              column(10,align="center",DT::dataTableOutput('dt')),
              column(1)
            ),
            fluidRow(
              h3("SNPs and Indels",align = "center"),
              column(12,align="center",plotOutput('snptype', height = "600px"))
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
