#! /usr/bin/env Rscript
#
#-------------------------------------------------------------------#
# Insert metabolome data from readXLS.R into MongoDB database
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

library(mongodb)
library(assertthat)

# Connect to MongoDB ---------------------------------------------------

mongo <- mongo.create()

assert_that(mongo.is.connected(mongo))
