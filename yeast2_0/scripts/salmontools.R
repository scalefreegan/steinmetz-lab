#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# Tools for tRNA-seq analysis with salmon
#
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2016"
.credits = "Aaron Brooks"
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"

library(wasabi)
library(sleuth)
library(readr)
library(rhdf5)

processSalmonH5 <- function(f) {
    dh5 = h5read(as.vector(f),"/")
    o = data.frame(id = dh5$aux$ids, est = dh5$est_counts, do.call(cbind,dh5$bootstrap),
      stringsAsFactors = F)
    H5close()
    return(o)
}

processEQ <- function(f) {
    fi = read_delim(f,delim = "\n",col_names = "")
    tx = as.integer(fi[1,1])
    eq = fi[2,1]
    cat(paste("Processed file with", tx, "transcripts in", eq, "equivalence classes.\n"),sep="")
    txs = fi[2:tx+2,1]
    eqs = fi[tx+3:dim(fi)[1],1]
    eq2tx = eqs %>% rowwise() %>% do({
      ss = strsplit(as.character(.),"\t")[[1]]
      tx.per.eq = ss[1]
      frag.per.eq = ss[length(ss)]
      ids = ss[2:(length(ss)-1)]
      nids = paste(unlist(txs[as.integer(ids),1]), collapse = " ")
      data.frame(tx.per.eq, frag.per.eq, txs = nids, stringsAsFactors = F)
    })
    eq2tx$id = seq(1,dim(eqs)[1])
    # o$tx2eq = do.call(rbind,lapply(seq(1,length(txs[[1]])),function(x){
    #   tx = txs[[1]][x]
    #   eqs = paste(which(grepl(as.character(tx), o$eq2tx$txs)), collapse = " ")
    #   data.frame(tx = tx, eqs)
    # }))
    return(eq2tx)
}

tx2eq <- function(eq2tx) {
    # get txs
    txs = unique(unlist(sapply(unlist(eq2tx$txs),function(x){strsplit(x," ")[[1]]})))
    tx2eq = do.call(rbind,lapply(txs, function(x){
        neqs = length(which(grepl(x, eq2tx$txs)))
        eqs = paste(which(grepl(x, eq2tx$txs)), collapse = " ")
        data.frame(tx = x, eqs, neqs)
    }))
}
