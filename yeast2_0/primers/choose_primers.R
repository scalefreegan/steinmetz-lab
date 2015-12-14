#! /usr/bin/env Rscript
# designed to be saved as .Rprofile in working dir
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/metabolome/processMetabData_allstrains.R")
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
.plot = FALSE

# Import packages ---------------------------------------------------
library(plyr)
library(dplyr)
library(reshape2)

# Load files --------------------------------------------------------
