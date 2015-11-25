.local = FALSE
if (system("hostname",intern=T) == "mac-steinmetz55.embl.de") {
  print("yes")
  .local = TRUE
} else {
  print(system("hostname"))
}

# Import packages ---------------------------------------------------

library(shiny)
library(clustQTL)
library(dplyr)
library(reshape2)
library(DT)
library(parallel)
library(GenomicRanges)
library(ggplot2)
library(stringr)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("org.Sc.sgd.db")

library(RCurl)
library(httr)
set_config( config( ssl_verifypeer = 0L ) )
library(funqtl)

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
  DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"
} else {
  DDIR = "/g/steinmetz/brooks/genphen/metabolome/qtls"
}


f = file.path(DDIR,"mQTL.rda")
if (!file.exists(f)) {
  # load and process data
  devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/mQTL_explorer/processData.R")
} else {
  load(f)
}

#load("//g/steinmetz/brooks/genphen/metabolome/qtls/plot_tab.rda")
