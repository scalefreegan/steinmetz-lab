
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
library(GenomicRanges)
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

shinyServer(function(input, output, session) {

  options(DT.options = list(pageLength = 5))
  
  gbrowse_link = function(chr,start,end,flanking = c(-2000,2000)) {
    start = start - abs(flanking[1])
    if (start < 0) {
      start = 1
    }
    end = end + abs(flanking[2])
    val = paste("chr",chr,":",start,"..",end, sep = "")
    sprintf('<a href="http://browse.yeastgenome.org/fgb2/gbrowse/scgenome/?name=%s" target="_blank" class="btn btn-primary">Go to gene</a>',
            val)
  }
  
  gbrowse_add = function(chr,start,end) {
    val = paste("chr",chr,":",start,"..",end, sep = "")
    sprintf('<a href="http://browse.yeastgenome.org/fgb2/gbrowse/scgenome/?add=%s" target="_blank" class="btn btn-danger">Add data</a>',
            val)
  }
  
  df = reactive({
    qtl_df = as.data.frame(do.call(rbind,lapply(names(clust_qtls),function(i){
      i_2 = gsub("/","%2",i)
      strand = GenomicRanges::strand(tx_3utr[i_2])
      chr = as.character(as.roman(as.numeric(gsub("chr","",GenomicRanges::seqnames(tx_3utr[i_2])))))
      start = IRanges::start(IRanges::ranges(tx_3utr[i_2]))
      end = IRanges::end(IRanges::ranges(tx_3utr[i_2]))
      type = GenomicRanges::elementMetadata(tx_3utr[i_2])$type
      i_peaks = findQTLPeaks(clust_qtls[[i]]$qtl, mrk, pcutoff = 10^-input$alpha)
      if (length(i_peaks) == 0 ) {
        data.frame(Gene = i, QTLs = length(i_peaks), Top_QTL = NA, 
                   Strand = strand, Chr = chr, Start = start, End = end, Type = type, 
                   GBrowse = paste(gbrowse_link(chr,start,end),gbrowse_add(chr,start,end)))
      } else {
        data.frame(Gene = i, QTLs = length(i_peaks), Top_QTL = max(i_peaks$p), 
                   Strand = strand, Chr = chr, Start = start, End = end, Type = type, 
                   GBrowse = paste(gbrowse_link(chr,start,end),gbrowse_add(chr,start,end)))
      }
    })))
    qtl_df[order(qtl_df$QTLs,qtl_df$Top_QTL,decreasing=T),]
  })
  
  # DATA TABLE
  output$dt = DT::renderDataTable(
    df(), server = TRUE, selection = "single",
    rownames = FALSE, extensions = 'Responsive', escape = FALSE)

  # REACTIVE VALUES
  values = reactiveValues(
    old_selection = NULL,
    marker = NULL,
    x = NULL
  )
  
  # MONITOR OLD SELECTION
  session$onFlush(once=FALSE, function(){
    isolate({ values$old_selection <- input$dt_rows_selected })
  })
  
  # OBSERVE CLICK
  observeEvent(input$plot_click, {
    marker = nearPoints(manhattan_data(),input$plot_click, xvar = "pos", yvar = "logp")
    if (dim(marker)[1] > 0) {
      marker = marker[which.max(marker$logp),"SNP"]
      values$marker = levels(marker)[marker]
    }
  })
  
  # OBSERVE BRUSH
  observeEvent(input$plot_brush, {
    marker = brushedPoints(manhattan_data(), input$plot_brush, xvar = "pos", yvar = "logp")
    if (dim(marker)[1] > 0) {
      marker = marker[which.max(marker$logp),"SNP"]
      values$marker = levels(marker)[marker]
    }
  })
  
  # OBSERVE DOUBLECLICK
  observeEvent(input$plot_dblclick, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      x <- brushedPoints(manhattan_data(), input$plot_brush, xvar = "pos", yvar = "logp")[,"SNP"]
      values$x <- levels(x)[x]
    } else {
      values$x <- NULL
    }
  })
  
  gene = reactive({
    s = input$dt_rows_selected
    if (length(s)) {
     levels(df()[s, "Gene"])[df()[s, "Gene"]]
    } else {
      NULL
    }
  })
  
  gene2 = reactive({
    s = input$dt_rows_selected
    if (length(s)) {
      gsub("/", "%2", df()[s, "Gene"])
    } else {
      NULL
    }
  })
  
  manhattan_data = reactive({
    if (!is.null(gene())) {
      qtls = clust_qtls[[gene()]]$qtl
      rownames(qtls) = names(mrk)
      # restrict based on 
      if (!is.null(values$x)) {
        qtls = qtls[values$x,]
        tmp_mrk = mrk[values$x]
      } else {
        tmp_mrk = mrk
      }
      clustQTL::plotManhattan(qtls, tmp_mrk, gene = gene2(), trx_annot = tx_3utr,
                              cutoff = 10^-input$alpha, qqman = TRUE, show = FALSE, 
                              suggestiveline = input$alpha, genomewideline = FALSE)
    }
  })
  
  output$manhattan = renderPlot({
    if (!is.null(gene())) {
      qtls = clust_qtls[[gene()]]$qtl
      rownames(qtls) = names(mrk)
      # restrict based on 
      if (!is.null(values$x)) {
        qtls = qtls[values$x,]
        tmp_mrk = mrk[values$x]
      } else {
        tmp_mrk = mrk
      }
      clustQTL::plotManhattan(qtls, tmp_mrk, gene = gene2(), trx_annot = tx_3utr,
                      cutoff = 10^-alpha, qqman = TRUE,
                      suggestiveline = input$alpha, genomewideline = FALSE)
    }
  })
  
  output$pqtl = renderPlot({
    if (!is.null(gene())) {
      data <- clust_qtls[[gene()]]$data
      
      # find peaks
      peaks =  clustQTL::findQTLPeaks(clust_qtls[[gene()]]$qtl, mrk,
                           pcutoff = 10^-input$alpha)
      
      # draw the qtl peak profile
      if(input$dt_rows_selected != values$old_selection || is.null(values$marker)) {
        print("normal")
        print(names(peaks)[1])
        clustQTL::plotPeakProfile(data, geno,
                        names(peaks)[1], peak_sigma = 2,
                        peak_threshold = 1)
      } else {
          clustQTL::plotPeakProfile(data, geno,
                values$marker, peak_sigma = 2,
                peak_threshold = 1)
      }
    }
  })
})
