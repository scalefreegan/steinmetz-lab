#-------------------------------------------------------------------#
# Load yeast genetics stuff
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

# Import packages ---------------------------------------------------

library(dplyr)
library(reshape2)
library(parallel)
library(GenomicRanges)
library(ggplot2)
library(funqtl)
library(stringr)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("org.Sc.sgd.db")
library(RCurl)
library(httr)
set_config( config( ssl_verifypeer = 0L ) )
# Global variables ---------------------------------------------------

id2name = id2name(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

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
