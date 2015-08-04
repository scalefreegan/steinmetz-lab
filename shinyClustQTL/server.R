
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

if (!is.null(local)) {
  load("~/Desktop/tmpdata/3TFill/clust_qtl_small.rda")
  load("~/Desktop/tmpdata/3TFill/geno_mrk.RData")
  load("~/Desktop/tmpdata/3TFill/tx_3utr.rda")
} else {
  load("/g/steinmetz/brooks/3prime_Tfill/clust_qtl.rda")
  #load("/g/steinmetz/brooks/3prime_Tfill/clust_qtl_1000.rda")
  load( "/g/steinmetz/brooks/genphen/qtl_endometabolome_23042015/geno_mrk.RData" )
  load("/g/steinmetz/brooks/3prime_Tfill/tx_3utr.rda")
}
# load qt_file
#
#qtl_genes = sort(unlist(sapply(clust_qtls,function(i){min(as.numeric(i$qtl[,2]))})))
names(clust_qtls) = gsub("%2","/",names(clust_qtls))

# genotypes/markers from chenchen. only applies to clust_qtl data set
# will be replaced in future

# supplies mrk and geno
alpha = 1e-5
qtl_df = as.data.frame(do.call(rbind,lapply(names(clust_qtls),function(i){
  #alpha = 10^-input$alpha
  i_peaks = findQTLPeaks(clust_qtls[[i]]$qtl, mrk, pcutoff = alpha)
  if (length(i_peaks) == 0 ) {
    data.frame(gene = i,n_qtls = length(i_peaks), top_qtl = NA)
  } else {
    data.frame(gene = i,n_qtls = length(i_peaks), top_qtl = max(i_peaks$p))
  }
})))

df = qtl_df[order(qtl_df$n_qtls,qtl_df$top_qtl,decreasing=T),]

shinyServer(function(input, output, session) {

  options(DT.options = list(pageLength = 5))

  #input = list()
  #input$alpha = 5

  output$dt = DT::renderDataTable(
    df, server = TRUE, selection = "single")

  selection = reactiveValues(
    old = "start"
  )
  
  session$onFlush(once=FALSE, function(){
    isolate({ selection$old<-input$dt_rows_selected })
  })
  
  gene = reactive({
    s = input$dt_rows_selected
    #input$plot_click = NULL
    #input$plot_brush = NULL
    print(s)
    print(selection$old)
    if (length(s)) {
     levels(df[s, "gene"])[df[s, "gene"]]
    } else {
      NULL
    }
  })
  
  gene2 = reactive({
    s = input$dt_rows_selected
    if (length(s)) {
      gsub("/", "%2", df[s, "gene"])
    } else {
      NULL
    }
  })
  
  manhattan_data = reactive({
    if (!is.null(gene())) {
      clustQTL::plotManhattan(clust_qtls[[gene()]]$qtl, mrk, gene = gene2(), trx_annot = tx_3utr,
                              cutoff = alpha, qqman = TRUE, show = FALSE, 
                              suggestiveline = -log10(alpha), genomewideline = FALSE)
    }
  })
  
  output$manhattan = renderPlot({
    if (!is.null(gene())) {
      qtls = clust_qtls[[gene()]]$qtl
      clustQTL::plotManhattan(qtls, mrk, gene = gene2(), trx_annot = tx_3utr,
                      cutoff = alpha, qqman = TRUE,
                      suggestiveline = -log10(alpha), genomewideline = FALSE)
    }
  })
  
  output$pqtl = renderPlot({
    if (!is.null(gene())) {
      data <- clust_qtls[[gene()]]$data
      
      # find peaks
      peaks =  clustQTL::findQTLPeaks(clust_qtls[[gene()]]$qtl, mrk,
                           pcutoff = alpha)
      
      # draw the qtl peak profile
      if( (is.null(input$plot_click) && is.null(input$plot_brush) ) || 
          input$dt_rows_selected != selection$old ) {
        print("normal")
        print(names(peaks)[1])
        clustQTL::plotPeakProfile(data, geno,
                        names(peaks)[1], peak_sigma = 2,
                        peak_threshold = 1)
      } else {
        if(is.null(input$plot_brush)) {
          print("click")
          marker = nearPoints(manhattan_data(),input$plot_click, xvar = "pos", yvar = "logp")
        } else {
          print("brush")
          marker = brushedPoints(manhattan_data(), input$plot_brush, xvar = "pos", yvar = "logp")
          }
        if (dim(marker)[1] > 0) {
          marker = marker[which.max(marker$logp),"SNP"]
          marker = levels(marker)[marker]
          print(marker)
          clustQTL::plotPeakProfile(data, geno,
                          marker, peak_sigma = 2,
                          peak_threshold = 1)
        }
      }
    }
  })
})
