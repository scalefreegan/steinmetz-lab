
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)
library(clustQTL)
library(dplyr)
library(DT)
library(parallel)
options("mc.cores"=4)
# load qt_file
#load("/g/steinmetz/brooks/3prime_Tfill/clust_qtl.rda")
load("~/Desktop/tmpdata/3TFill/clust_qtl_small.rda")
#qtl_genes = sort(unlist(sapply(clust_qtls,function(i){min(as.numeric(i$qtl[,2]))})))
names(clust_qtls) = gsub("%2","/",names(clust_qtls))

# genotypes/markers from chenchen. only applies to clust_qtl data set
# will be replaced in future
#load( "/g/steinmetz/brooks/genphen/qtl_endometabolome_23042015/geno_mrk.RData" )
load("~/Desktop/tmpdata/3TFill/geno_mrk.RData")
load("~/Desktop/tmpdata/3TFill/tx_3utr.rda")
# supplies mrk and geno
alpha = 2
qtl_df = as.data.frame(do.call(rbind,mclapply(names(clust_qtls),function(i){
  #alpha = 10^-input$alpha
  i_peaks = findQTLPeaks(clust_qtls[[i]]$qtl, mrk, pcutoff = alpha)
  if (length(i_peaks) == 0 ) {
    data.frame(gene=i,n_qtls = length(i_peaks), top_qtl = NA)
  } else {
    data.frame(gene=i,n_qtls = length(i_peaks), top_qtl = max(i_peaks$p))
  }
})))

df = qtl_df[order(qtl_df$n_qtls,qtl_df$top_qtl,decreasing=T),]

shinyServer(function(input, output, session) {

  options(DT.options = list(pageLength = 5))

  #input = list()
  #input$alpha = 5

  output$dt = DT::renderDataTable(
    df, server = TRUE, selection = "single")

  output$manhattan = renderPlot({
    s = input$dt_rows_selected
    if (length(s)) {
      gene =  df[s, "gene"]
      gene2 =  gsub("/", "%2", df[s, "gene"])
      p = plotManhattan(clust_qtls[[gene]]$qtl, mrk, gene = gene2, trx_annot = tx_3utr,
                              cutoff = alpha)
      p
    }
  })
    
  output$pqtl = renderPlot({
    s = input$dt_rows_selected
    if (length(s)) {
      gene = df[s,"gene"]
      print(gene)
      data <- clust_qtls[[gene]]$data
      
      # find peaks
      peaks = findQTLPeaks(clust_qtls[[gene]]$qtl, mrk,
                           pcutoff = alpha)
      
      # draw the qtl peak profile
      plotPeakProfile(data, geno,
                      names(peaks)[1], peak_sigma = 2,
                      peak_threshold = 1)
    }
  })
  
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    #mrk2 = as.data.frame(format4manhattan( clust_qtls[[gene]]$qtl,mrk ))
    paste0(
      "click: ", xy_str(input$plot_click),
      "brush: ", xy_range_str(input$plot_brush)
    )
  })
  
})
