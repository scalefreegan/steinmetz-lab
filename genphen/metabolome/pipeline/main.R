#! /usr/bin/env Rscript
#
#-------------------------------------------------------------------#
# Run the entire pipeline
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

library(devtools)

# Definitions ---------------------------------------------------
DATADIR = "/g/steinmetz/project/GenPhen/data/"
DBNAME = "genphen"

# Load components ---------------------------------------------------

source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/metabolome/pipeline/readXLS.R")
source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/metabolome/pipeline/insertMongoDB.R")

# Process data ---------------------------------------------------

f1 = paste(DATADIR, "endometabolome/data/Endometabolome_1B_46B_sorted by cultivation phase.xlsx", sep="")
f2 = paste(DATADIR, "endometabolome/data/Endometabolome_46B_sorted by cultivation time.xlsx", sep="")

thisdata = processData(f1, f2)
