#! /usr/bin/env Rscript
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/integration/integration.R")

#-------------------------------------------------------------------#
# GenPhen segregant integrative Analysis
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

# library( plotly )
# set_credentials_file( "scalefreegan", "5xyvvaleim" )
library( qtl )
library( GenomicRanges )
library( rtracklayer )
library( plyr )
library( parallel )
options("mc.cores"=24)
library( ggplot2 )
library( ggbio )
library( LSD )
library( DESeq2 )
library(dplyr)
library(reshape2)

# source rQTL utilities
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/QTL/rQTL.R")

.plot = FALSE

#-------------------------------------------------------------------#
# Load data
#-------------------------------------------------------------------#

DDIR = "/g/steinmetz/project/GenPhen/data/endometabolome"
MDIR = "/g/steinmetz/brooks/genphen/metabolome"
EDIR = "/g/steinmetz/brooks/genphen/transcriptome/qtl"

load(file.path(DDIR,"data/endometabolite_full_12102015.rda"))
load(file.path("/g/steinmetz/brooks/yeast/genomes/S288CxYJM789/genotypes_S288c_R64.rda"))
load(file.path(MDIR,"qtls/mQTLs_comball_funqtl_2014.rda"))
load(file.path(EDIR,"eQTL.rda"))

#-------------------------------------------------------------------#
# mQTL/eQTL correlation
#-------------------------------------------------------------------#
