# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#
#-------------------------------------------------------------------#
# Shiny interface for ploting/exploring mQTLs
# 
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "WTFPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.plot = FALSE

# Import packages ---------------------------------------------------

library(shiny)
library(clustQTL)
library(dplyr)
library(reshape2)
library(DT)
library(parallel)
library(GenomicRanges)
library(ggplot2)
library("org.Sc.sgd.db")
orf2name = org.Sc.sgdGENENAME

addResourcePath('data', "/Users/brooks/Sites/JBrowse-1.11.6/data")

# tmp resource location / will be changed
DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"

# composite rda
f = file.path(DDIR,"mQTL.rda")
if (!file.exists(f)) {
  # load and process data

} else {
  load(f)
}


# try to write

if (!is.null(local)) {
  load("~/Desktop/tmpdata/3TFill/clust_qtl_small.rda")
  load("~/Desktop/tmpdata/3TFill/geno_mrk.RData")
  load("~/Desktop/tmpdata/3TFill/tx_3utr.rda")
  #load("~/Desktop/tmpdata/3TFill/mergedCounts.rda")
} else {
  #load("/g/steinmetz/brooks/3prime_Tfill/clust_qtl.rda")
  load("/g/steinmetz/brooks/3prime_Tfill/clust_qtl_1000.rda")
  load( "/g/steinmetz/brooks/genphen/qtl_endometabolome_23042015/geno_mrk.RData" )
  load("/g/steinmetz/brooks/3prime_Tfill/tx_3utr.rda")
  #load("/g/steinmetz/project/GenPhen/data/3tagseq/all/mergedCounts.rda")
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
    val = paste("chr",as.roman(chr),":",start,"..",end, sep = "")
    #sprintf('<a href="http://browse.yeastgenome.org/fgb2/gbrowse/scgenome/?name=%s" target="_blank" class="btn btn-primary">Go to gene</a>',
    #        val)
    sprintf('<a href="http://localhost/~brooks/JBrowse-1.11.6/?loc=%s" target="_blank" class="btn btn-primary">Go to gene</a>',
            val)
  }
  
  gbrowse_heattrack = function(data) {
    
    STORETEMPLATE = "{'%s':{'type':'JBrowse/Store/SeqFeature/BigWig','urlTemplate':'%s'}}"
    
    TRACKTEMPLATE = "{'label':'%s','key':'%s','urlTemplate':'%s','type':'JBrowse/View/Track/Wiggle/Density','bicolor_pivot':'mean','style':{'pos_color':'purple','neg_color':'green'}}"
    
    TRACKTEMPLATE = "{'label':'%s','key':'%s','urlTemplate':'%s','type':'JBrowse/View/Track/Wiggle/Density','storeClass':'JBrowse/Store/SeqFeature/BigWig','bicolor_pivot':'mean','style':{'pos_color':'purple','neg_color':'green'}}"
    
    
    
    this_store = paste("addStores=",
                       URLencode(paste(sprintf(STORETEMPLATE,"url","t1.bw"),sep=""), reserved = T), sep = "")
    this_tracks = paste("addTracks=",
                      URLencode(paste("[",sprintf(TRACKTEMPLATE,"iso1","S288c_Isoforms","url"),"]",sep=""),reserved = T),
                      sep = "")
    this_tracks = paste("addTracks=",
                        URLencode(paste("[",sprintf(TRACKTEMPLATE,"iso1","S288c_Isoforms","t1.bw"),"]",sep=""),reserved = T),
                        sep = "")
    o = paste(this_store,"&",this_tracks,sep="")     
    #o = paste(this_tracks,sep="")     
  }
  
  
  
  makeBigWig = function(data, this_mrk, geno) {
    d1 = colSums(data[intersect(names(which(geno[this_mrk,]==1)),rownames(data)),])
    d2 = colSums(data[intersect(names(which(geno[this_mrk,]==2)),rownames(data)),])
    chr = df()[s,"`Chr"]
    seqinfo = Seqinfo(genome="sacCer3")
    g1 = GenomicRanges::GRanges(
      #seqnames = Rle(seqnames(chr_info[chr]),seqlengths(chr_info[chr])), 
      seqnames = Rle(seqnames(chr_info[chr]),1),
      ranges = IRanges(as.numeric(names(d1)), end = as.numeric(names(d1))),
      strand = "+",
      score = d1,
      seqinfo = seqinfo
      )
    #seqlengths(g1) = seqlengths(chr_info[chr])
    #seqlengths(g1) = length(d1)
    rtracklayer::export.bw(object = g1, con = "/Users/brooks/Sites/JBrowse-1.11.6/data/t1.bw", compress=T)
    g2 = GenomicRanges::GRanges(
      #seqnames = Rle(seqnames(chr_info[chr]),seqlengths(chr_info[chr])), 
      seqnames = Rle(seqnames(chr_info[chr]),1),
      ranges = IRanges(as.numeric(names(d2)), end = as.numeric(names(d2))),
      strand = "+",
      score = d2,
      seqinfo = seqinfo
    )
    #seqlengths(g2) = seqlengths(chr_info[chr])
    rtracklayer::export.bw(object = g2, con = "/Users/brooks/Sites/JBrowse-1.11.6/data/t2.bw", compress=T)
  }
     
  df = reactive({
    qtl_df = as.data.frame(do.call(rbind,lapply(names(clust_qtls),function(i){
      i_2 = gsub("/","%2",i)
      i_out = orf2name[[i]]
      if (is.null(i_out) || is.na(i_out)) {
        i_out = i
      }
      strand = GenomicRanges::strand(tx_3utr[i_2])
      chr = as.character(as.numeric(gsub("chr","",GenomicRanges::seqnames(tx_3utr[i_2]))))
      start = IRanges::start(IRanges::ranges(tx_3utr[i_2]))
      end = IRanges::end(IRanges::ranges(tx_3utr[i_2]))
      type = GenomicRanges::elementMetadata(tx_3utr[i_2])$type
      i_peaks = findQTLPeaks(clust_qtls[[i]]$qtl, mrk, pcutoff = 10^-input$alpha)
      if (length(i_peaks) == 0 ) {
#         data.frame(Name = i_out, ORF = i, QTLs = length(i_peaks), Top_QTL = NA, 
#                    Strand = strand, Chr = chr, Start = start, End = end, Type = type, 
#                    JBrowse = paste(gbrowse_link(chr,start,end)))
          data.frame(Name = i_out, ORF = i, QTLs = length(i_peaks), Top_QTL = NA, 
                 Strand = strand, Chr = chr, Start = start, End = end, Type = type)
      } else {
#         data.frame(Name = i_out, ORF = i, QTLs = length(i_peaks), Top_QTL = max(i_peaks$p), 
#                    Strand = strand, Chr = chr, Start = start, End = end, Type = type, 
#                    JBrowse = paste(gbrowse_link(chr,start,end)))
          data.frame(Name = i_out, ORF = i, QTLs = length(i_peaks),
              Top_QTL = max(i_peaks$p), Strand = strand, Chr = chr, Start =
              start, End = end, Type = type)
      }
    })))
    qtl_df[order(qtl_df$QTLs,qtl_df$Top_QTL,decreasing=T),]
  })
  
  # DATA TABLE
  output$dt = DT::renderDataTable(
    df(), server = TRUE, selection = "single",
    rownames = FALSE, extensions = 'Responsive', escape = FALSE)

  # STATIC VALUES
  out$metabolites = 
  
  # REACTIVE VALUES
  values = reactiveValues(
    old_selection = NULL,
    marker = NULL,
    x = NULL,
    link = NULL
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
     levels(df()[s, "ORF"])[df()[s, "ORF"]]
    } else {
      NULL
    }
  })
  
  gene2 = reactive({
    s = input$dt_rows_selected
    if (length(s)) {
      gsub("/", "%2", df()[s, "ORF"])
    } else {
      NULL
    }
  })
  
  #output$link = renderPrint("hi")
  output$link = renderText({
    s = input$dt_rows_selected
    print(s)
    if (length(s)) {
      #print(df()[s,])
      chr = levels(df()[s, "Chr"])[df()[s, "Chr"]]
      start = df()[s, "Start"]-1000
      end = df()[s, "End"]+1000
      val = paste("chr",as.roman(chr),"%3A",start,"..",end, sep = "")
      val = paste('http://localhost/~brooks/JBrowse-1.11.6/?loc=',val,sep="")
    } else {
      val = "http://localhost/~brooks/JBrowse-1.11.6/"
    }
    paste("<div style='width: 100%; height: 600px'><iframe style='border: 1px solid black' src='", val ,"'width='100%' height='100%'></iframe></div>",sep="")
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
