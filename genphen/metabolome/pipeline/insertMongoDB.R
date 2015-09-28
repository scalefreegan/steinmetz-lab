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

library(rmongodb)
library(assertthat)

# Connect to MongoDB ---------------------------------------------------

mongoConnect = function(host = NULL, db = NULL) {
  if (is.null(host) || is.null(db)) {
    cat("Must supply host (host=) and database (db=) to connect to MongoDB\n")
    return(NULL)
  }
  mongo <- mongo.create(host = , db = db)
  if (assert_that(mongo.is.connected(mongo))) {
    return(mongo)
  } else {
    cat("ERROR: Could not connect to MongoDB\n")
    return(NULL)
  }


}
