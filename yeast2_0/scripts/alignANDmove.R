#!/usr/bin/env Rscript

#-------------------------------------------------------------------#
# Align query and reference sequence. Move the origin of
# query sequence to match origin of reference
# Written to align PacBio assemblies w/ short read assemblies
#
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2016"
.credits = "Aaron Brooks"
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"

# Import packages ---------------------------------------------------
if (system('hostname',intern=T) == "taiji.embl.de") {
    .libPaths("~/R/x86_64-redhat-linux-gnu-library/3.2/")
}

library(optparse)

option_list = list(
  make_option(c("-r", "--ref"), type="character", default=NULL,
              help="Reference sequence. Fasta format.", metavar="character"),
	make_option(c("-q", "--query"), type="character",
              help="Query sequence. Fasta format.", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="ref_query",
              help="Output name for MUMMER alignment files. [default= %default]",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$ref)) {
  print_help(opt_parser)
  stop("Reference sequence in fasta format required.n", call.=FALSE)
} else if (is.null(opt$query)) {
  print_help(opt_parser)
  stop("Query sequence in fasta format required.n", call.=FALSE)
}

mummer_command = paste("nucmer --prefix=", opt$prefix, "--coords", opt$ref, opt$query, sep = " ")

system(mummer_command, wait = TRUE)
