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
library(stringr)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("org.Sc.sgd.db")

# Global variables ---------------------------------------------------

id2name = id2name(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
type = "mlod"

# gene name map
x = org.Sc.sgdGENENAME
keys <- mappedkeys(x)
# Convert to a list
gname <- as.list(x[keys])

# short description map
x = org.Sc.sgdALIAS
keys <- mappedkeys(x)
# Convert to a list
dname <- as.list(x[keys])

# short description map
x = org.Sc.sgdDESCRIPTION
keys <- mappedkeys(x)
# Convert to a list
dname_long <- as.list(x[keys])

# Web resources ---------------------------------------------------

addResourcePath('data', "/Users/brooks/Sites/JBrowse-1.11.6_mQTL/data")

# Misc material ---------------------------------------------------

# composite rda
DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"

f = file.path(DDIR,"mQTL.rda")
if (!file.exists(f)) {
  # load and process data
  devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/mQTL_explorer/processData.R")
} else {
  load(f)
}

shinyServer(function(input, output, session) {

  options(DT.options = list(pageLength = 10, searching = TRUE))
  
  gbrowse_link = function(chr,start,end,flanking = c(-2000,2000)) {
    start = start - abs(flanking[1])
    if (start < 0) {
      start = 1
    }
    end = end + abs(flanking[2])
    val = paste("chr",as.roman(chr),":",start,"..",end, sep = "")
    #sprintf('<a href="http://browse.yeastgenome.org/fgb2/gbrowse/scgenome/?name=%s" target="_blank" class="btn btn-primary">Go to gene</a>',
    #        val)
    sprintf('<a href="http://localhost/~brooks/JBrowse-1.11.6_mQTL/?loc=%s" target="_blank" class="btn btn-primary">Go to gene</a>',
            val)
  }
  
  df = reactive({
    
    # select chromosomes
    chrs = unique(data[[input$m]]$qtl[data[[input$m]]$qtl[,type]>=summary(data[[input$m]]$permout[,type],input$co/100)[1],"chr"]) 
    chrs = levels(chrs)[chrs]
    
    lodcolumn = if(type=="mlod"){ 2 } else { 1 }
    qtl_intervals = list()
    if (length(chrs)>0) {
      for (i in chrs) {
        qtl_intervals[[i]] = try(mrk[rownames(bayesint(data[[input$m]]$qtl, chr = str_pad(i, 2, pad = "0"), prob=input$bci/100, lodcolumn=lodcolumn))],silent = T)
        if (class(qtl_intervals[[i]])=="try-error") {
          qtl_intervals[[i]] = NULL
        } else {
          nn = sapply(as.character(seqnames(qtl_intervals[[i]])),function(i){
            paste(substr(i,1,3),as.roman(substr(i,4,5)),sep="")
          })
          qtl_intervals[[i]] = renameSeqlevels(qtl_intervals[[i]],nn)
          qtl_intervals[[i]] = keepSeqlevels(qtl_intervals[[i]],unique(nn))
          qtl_intervals[[i]] = range(qtl_intervals[[i]])
          qtl_intervals[[i]] = as.data.frame(cdsByOverlaps(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,qtl_intervals[[i]], type = "any", columns = "gene_id"))
        }
      }
    }
    qtl_df = do.call(rbind,qtl_intervals)
    if (length(qtl_df) != 0) {
      qtl_df$gene_id = unlist(qtl_df$gene_id)
      gname_t = unlist(gname[unlist(qtl_df$gene_id)])
      gname_t = data.frame(gene_id = names(gname_t), name = gname_t)
      dname_t = unlist(dname[unlist(qtl_df$gene_id)])
      dname_t = data.frame(gene_id = names(dname_t), alias = dname_t)
      dname_t_long = unlist(dname_long[unlist(qtl_df$gene_id)])
      dname_t_long = data.frame(gene_id = names(dname_t_long), desc = dname_t_long)
      qtl_df = merge(qtl_df,gname_t,by="gene_id",sort=F,all.x=T)
      qtl_df = merge(qtl_df,dname_t,by="gene_id",sort=F,all.x=T)
      qtl_df = merge(qtl_df,dname_t_long,by="gene_id",sort=F,all.x=T)
      qtl_df = qtl_df[,c("gene_id","name","seqnames","start","end","strand","alias","desc")]
      colnames(qtl_df) = c("Sys.Name","Name","Chr","Start","End","Strand","Alias","Desc")
      qtl_df = qtl_df[!duplicated(qtl_df),]
    }
  })
  
  # DATA TABLE
  output$dt = DT::renderDataTable(
    df(), server = TRUE, selection = "single", 
    rownames = FALSE, extensions = 'Responsive', escape = FALSE)

  # REACTIVE VALUES
  values = reactiveValues(
    old_selection = NULL,
    link = NULL
  )
  
  # MONITOR OLD SELECTION
  session$onFlush(once=FALSE, function(){
    isolate({ values$old_selection <- input$dt_rows_selected })
  })
  
  #output$link = renderPrint("hi")
  output$link = renderText({
    s = input$dt_rows_selected
    if (length(s)) {
      #print(df()[s,])
      chr = levels(df()[s, "Chr"])[df()[s, "Chr"]]
      start = df()[s, "Start"]-1000
      end = df()[s, "End"]+1000
      val = paste(chr,"%3A",start,"..",end, sep = "")
      val = paste('http://localhost/~brooks/JBrowse-1.11.6_mQTL/?loc=',val,sep="")
    } else {
      val = "http://localhost/~brooks/JBrowse-1.11.6_mQTL/"
    }
    paste("<div style='width: 100%; height: 600px'><iframe style='border: 1px solid black' src='", val ,"'width='100%' height='100%'></iframe></div>",sep="")
  })
  
  reformatQTL = reactive({
    qtls = data[[input$m]]$qtl
    names(qtls)[names(qtls) == type] = "pval"
    qtls[,"pval"]  = 10^-qtls[,"pval"]
    qtls
  })
  
  alpha_5 = reactive({
    summary(data[[input$m]]$permout[,type],.05)
  })
  
  alpha_10 = reactive({
    summary(data[[input$m]]$permout[,type],input$co/100)
  })
  
  ymax = reactive({
    max(max(alpha_5(), alpha_10()),max(-log10(reformatQTL()[,"pval"])))
  })
  
  manhattan_data = reactive({
    clustQTL::plotManhattan(reformatQTL(), mrk, qqman = TRUE, show = FALSE,
                            suggestiveline = alpha_5(), genomewideline = alpha_10(), ylab = "LOD", ylim = c(0,ymax()+2))
  })
  
  output$manhattan = renderPlot({
      clustQTL::plotManhattan(reformatQTL(), mrk, qqman = TRUE, show = TRUE,
                            suggestiveline = alpha_5(), genomewideline = alpha_10(), ylab = "LOD", ylim = c(0,ymax()+2))
      legend("topright", y.leg[i], c("5% FDR",paste(rownames(alpha_10())[1],"FDR")), lty = c(1, 1), col = c("blue", "red"))
  })
  
})
