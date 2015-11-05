# ! /usr/bin/env Rscript
# designed to be saved as .Rprofile in working dir
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/growth/growthQTL.R")
#-------------------------------------------------------------------#
# Detect growthQTLs for GenPhen Data from
# growth data from Nicole (metabolomics experiments)
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
library(ggplot2)
#library(plotly)
library(plyr)
library(dplyr)
library(reshape2)
library(LSD)
library(qtl)
library(pheatmap)
library(funqtl)
library(parallel)
options(mc.cores = 24)
library(snow)
#library(clustQTL)
#devtools::install_github("scalefreegan/steinmetz-lab/clustQTL")

# source rQTL utilities
devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/QTL/rQTL.R")

endo_f = "/g/steinmetz/project/GenPhen/data/endometabolome/data/endometabolite_full_12102015.rda"

# load complete endometabolite file
load(endo_f)
