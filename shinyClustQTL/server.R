
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
library(clustQTL)

load("/g/steinmetz/brooks/3prime_Tfill/clust_qtl.rda")
qtl_genes = sort(unlist(sapply(clust_qtls,function(i){min(as.numeric(i$qtl[,2]))})))

shinyServer(function(input, output) {
   
  output$qtlPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    data <- clust_qtl[[input$gene]]$data

    # find peaks
    peaks = findQTLPeaks(clust_qtl[[input$gene]]$qtl, markers_yeast, 
                         pcutoff = .01, peak_sigma = 25, peak_threshold=1)
    
    # draw the qtl peak profile
    plotPeakProfile(data, genotypes_yeast, 
                    names(peaks)[1], peak_sigma = 2, 
                    peak_threshold = 1)
    
  })
  
})
