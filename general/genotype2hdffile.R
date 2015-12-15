#! /usr/bin/env Rscript
#
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/metabolome/processMetabData_allstrains.R")
#-------------------------------------------------------------------#
# Process Metabolome Data Nicole on 07.07.2015
# ! Transform into a computable format !
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.plot = FALSE

# Import packages ---------------------------------------------------

library(dplyr)
library(reshape2)
library(parallel)
library(GenomicRanges)
library(ggplot2)
library(stringr)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("org.Sc.sgd.db")
library(rtracklayer)
library(jsonlite)
library(RCurl)
library(httr)
library(VariantAnnotation)

VDIR = "/Users/brooks/Documents/steinmetz_local/yeast/genomes/S288CxYJM789"
load(file.path(VDIR,"yjm789snpsANDindels_info.rda"))
load("/Users/brooks/Documents/git/steinmetz-lab/general/genotypes_S288c_R64.rda")

# Write geno file ---------------------------------------------------

# Try to match sequences in var_info and mrk - make table rownames = combined marker name & variant, colnames = strains
options(scipen=999)
S288CxYJM789 = geno
varNames = rownames(var_info)
varNames = data.frame(do.call(rbind,lapply(varNames, function(i){
  o = unlist(sapply(strsplit(i,":")[[1]],function(j){strsplit(j,"_")[[1]]}))
  names(o) = c("chr","pos","var")
  return(o)
  })),stringsAsFactors=F)
varNames$pos = as.numeric(varNames$pos)
gnames = unlist(mclapply(seq(1,dim(mrk)[1]),function(i){
  #print(i)
  entry = mrk[i,]
  match = filter(varNames,pos<=entry$end,pos>=entry$start)
  if (dim(match)[1]>0) {
    # keep only one match
    match = match[1,]
    var1 = strsplit(match$var,split="/")[[1]][1]
    var2 = strsplit(match$var,split="/")[[1]][2]
  } else {
    var1 = "N"
    var2 = "N"
  }
  return(paste(entry$start,entry$chr,entry$end,var1,var2,sep="_"))
}))

geno2write = geno
rownames(geno2write) = gnames
write.table(geno2write, row.names = T, col.names = T,
  file = "/Users/brooks/Documents/git/steinmetz-lab/general/S288CxYJM789.txt", quote = F)

# Write pheno file ---------------------------------------------------
endo_f = "/g/steinmetz/project/GenPhen/data/endometabolome/data/endometabolite_full_12102015.rda"
load(endo_f)
endodata = filter(endometabolite, time_format=="relative")
m = acast(endodata,formula=strain~metabolite+time,value.var="value.log2",fun.aggregate=mean)
rownames(m) = str_pad(rownames(m),width=3,side="left",pad="0")

write.table(m, row.names = T, col.names = T,
  file = "/Users/brooks/Documents/git/steinmetz-lab/general/phenotype_metabolites.txt", quote = F)
