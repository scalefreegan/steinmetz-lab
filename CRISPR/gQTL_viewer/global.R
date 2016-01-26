.local = FALSE
if (system("hostname",intern=T) == "mac-steinmetz55.embl.de" || system("hostname",intern=T) == "interzone.local") {
  print("yes")
  .local = TRUE
} else {
  print(system("hostname"))
}

# Import packages ---------------------------------------------------

library(shiny)
library(dplyr)
library(reshape2)
library(DT)
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
set_config( config( ssl_verifypeer = 0L ) )
library(VariantAnnotation)


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

#addResourcePath('data', "/var/www2/html/mQTL/data")

# Misc material ---------------------------------------------------

# composite rda
if (.local) {
  DDIR = "/Users/brooks/Documents/git/steinmetz_local/genphen/metabolome"
} else {
  DDIR = "/g/steinmetz/brooks/git/steinmetz-lab/CRISPR"
}


# tmp resource location / will be changed
if (.local) {
  DDIR = "/Users/brooks/Documents/git/steinmetz-lab/CRISPR"
  VDIR = "/Users/brooks/Documents/steinmetz_local/yeast/genomes/S288CxYJM789"
} else {
  DDIR = "/g/steinmetz/brooks/git/steinmetz-lab/CRISPR"
  VDIR = "/g/steinmetz/brooks/yeast/genomes/S288CxYJM789"
}

gQTL = read.delim(file.path(DDIR, "data/journal.pgen.1003803.s016.TXT"),sep="\t", skip=5)
load(file.path(VDIR,"yjm789snpsANDindels_info.rda"))

