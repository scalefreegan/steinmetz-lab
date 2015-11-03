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
.makeBW = FALSE

# Import packages ---------------------------------------------------
library(funqtl)
library(GenomicRanges)
library(rtracklayer)

# tmp resource location / will be changed
DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"

load(file.path(DDIR,"endometabolite_full_12102015.rda"))
load(file.path(DDIR,"genotypes_S288c_R64.rda"))
load(file.path(DDIR,"mQTLs_comball_funqtl_2014.rda"))

data = lapply(seq(1,length(mQTLs_funqtl_2014)), function(i){
  if (class(mQTLs_funqtl_2014[[i]])=="try-error") {
    return(NULL)
  } else {
    #print(i)
    o = list()
    o[["qtl"]] = mQTLs_funqtl_2014[[i]][["qtls_alt"]]
    o[["permout"]] = mQTLs_funqtl_2014[[i]]$permout
    return(o)
  }
})
names(data) = names(mQTLs_funqtl_2014)

# remove NULLs
data = data[!sapply(data,is.null)]

if (.makeBW) {
  nn = sapply(as.character(levels(seqnames(mrk))),function(i){
    paste(substr(i,1,3),as.roman(substr(i,4,5)),sep="")
  })
  mrk_mod = renameSeqlevels(mrk,nn)
  mrk_mod = keepSeqlevels(mrk_mod,unique(nn))
  mrk_mod$type = NULL
  mrk_mod$score = 1
  g = GenomicRanges::GRanges(
    seqnames = seqnames(mrk_mod),
    ranges = ranges(mrk_mod),
    strand = strand(mrk_mod),
    score = data[[1]]$qtl[,3],
    seqinfo = seqinfo
  )
  rtracklayer::export.bedGraph(object = g, con = "/Users/brooks/Sites/JBrowse-1.11.6_mQTL/data/t1.bw")
}
