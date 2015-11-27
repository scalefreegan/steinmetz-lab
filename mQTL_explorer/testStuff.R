devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/mQTL_explorer/global.R")

input = list()
input$m = "AKG"
input$co = 15
input$bci = 95
type = "mlod"

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