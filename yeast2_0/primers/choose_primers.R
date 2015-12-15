#! /usr/bin/env Rscript
# designed to be saved as .Rprofile in working dir
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/yeast2_0/primers/choose_primers.R")
#-------------------------------------------------------------------#
# Choose primers for strain validation
#
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.write = FALSE

# Import packages ---------------------------------------------------
library(plyr)
library(dplyr)
library(reshape2)

# Load files --------------------------------------------------------
GITHUBDIR = "http://scalefreegan.github.io/steinmetz-lab/yeast2_0/primers/"
synIII = read.delim(url(paste(GITHUBDIR, "synIII.tsv", sep = "")),header=T,sep="\t",stringsAsFactors=F)
synVI = read.delim(url(paste(GITHUBDIR, "synVI.tsv", sep = "")),header=T,sep="\t",stringsAsFactors=F)
synIXR = read.delim(url(paste(GITHUBDIR, "synIXR.tsv", sep = "")),header=T,sep="\t",stringsAsFactors=F)

# Process files --------------------------------------------------------
pSelect = function(data,n, customName = "SYNIII", cNames = c("X5..3..forward.Wild.Type","X5..3..reverse.Wild.Type","X5..3..forward.Synthetic","X5..3..reverse.Synthetic")) {
  # assume data order is linear along chr
  x = round(seq(1,dim(data)[1],length.out=n))
  s_data = data[x,]
  o = do.call(rbind,lapply(seq(1,dim(s_data)[1]),function(i){
    do.call(rbind,lapply(cNames,function(j){
      if (j=="X5..3..forward.Wild.Type") {
        d = "F"
        type = "WT"
      } else if (j=="X5..3..reverse.Wild.Type"){
        d = "R"
        type = "WT"
      } else if (j=="X5..3..forward.Synthetic") {
        d = "F"
        type = "SYN"
      } else if (j=="X5..3..reverse.Synthetic") {
        d = "R"
        type = "SYN"
      }
      nn = paste(customName,type,s_data[i,"ORF.Amp."],s_data[i,"Pair"],d,sep="_")
      return(cbind(nn,s_data[i,j]))
      }))
    }))
  return(o)
}

synIII_f = pSelect(synIII,5,customName="SYNIII")
synVI_f = pSelect(synVI,5,customName="SYNVI")
synIXR_f = pSelect(synIXR,5,customName="SYNIXR")
to_write = rbind(synIII_f,synVI_f,synIXR_f)
colnames(to_write) = c("Primer_Name","Sequence")

# Write files --------------------------------------------------------

if (.write) {
  write.table(to_write,file="~/Documents/git/steinmetz-lab/yeast2_0/primers/ordered.txt",quote=F,sep="\t",row.names=F)
}
