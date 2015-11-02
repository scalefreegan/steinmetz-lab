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
library(funqtl)
library("org.Sc.sgd.db")
orf2name = org.Sc.sgdGENENAME


addResourcePath('data', "/Users/brooks/Sites/JBrowse-1.11.6/data")

# tmp resource location / will be changed
DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"

# composite rda
f = file.path(DDIR,"mQTL.rda")
if (!file.exists(f)) {
  # load and process data
  devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/mQTL_explorer/processData.R")
} else {
  load(f)
}


# try to write

if (!is.null(local)) {
  load("~/Desktop/tmpdata/3TFill/clust_qtl_small.rda")
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
  
  output$manhattan = renderPlot({
      qtls = data[[input$m]]$qtl
      type = "mlod"
      names(qtls)[names(qtls) == type] = "pval"
      qtls[,"pval"]  = 10^-qtls[,"pval"]
      alpha_5 = data[[input$m]]$permout["5%",type]
      alpha_10 = data[[input$m]]$permout["10%",type]
      ymax = max(alpha_10,max(-log10(qtls[,"pval"])))
      clustQTL::plotManhattan(qtls, mrk, qqman = TRUE,
                      suggestiveline = alpha_5, genomewideline = alpha_10, ylab = "LOD", ylim = c(0,ymax+2))
  })
  
})
