
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)
library(markdown)

shinyUI(navbarPage("mQTL Explorer", selected = "Table Viewer",
         tabPanel("Table Viewer",id = "inspect", fluidPage(
            title = "mQTL Explorer",
            hr(),
            fluidRow(
              column(12,
                 fluidRow(
                   column(4,selectInput("m", " Metabolite", choices = names(data), selected = 1)),
                   column(4,sliderInput("co", " Custom FDR (%):,", max = 50, min = 0, value = 5)),
                   column(4,sliderInput("bci", " Bayesian Confidence Interval (%):,", max = 100, min = 0, value = 95))
                 )
                )
              ),
            fluidRow(
              column(12,
                 fluidRow(
                   column(12, htmlOutput('manhattan'))
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
