
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
            #h1("Poly(A) Isoform Viewer"),
            hr(),
            fluidRow(
              column(6,
                 fluidRow(
                   column(12,selectInput("m", " Metabolite", choices = names(data), selected = 1))
                 ),
                fluidRow(
                  column(6,sliderInput("bci", " Bayesian Confidence Interval (%):,", max = 100, min = 0, value = 50))
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
