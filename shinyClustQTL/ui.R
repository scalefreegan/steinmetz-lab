
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
qtl_genes = sort(unlist(sapply(clust_qtls,function(i){min(as.numeric(i$qtl[,2]))})))

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("3' UTR QTL Visualizer"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    selectInput(
      "gene", "Gene", choices = names(qtl_genes)
      )
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("qtlPlot")
  )
))
