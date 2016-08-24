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

# load  packages
library(ggplot2);
library(tidyr);
library(Biostrings);
library(seqinr);
library(optparse);

option_list = list(
  make_option(c("-r", "--ref"), type="character", default=NULL,
              help="Reference sequence. Fasta format.", metavar="character"),
	make_option(c("-q", "--query"), type="character",
              help="Query sequence. Fasta format.", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Name of re-aligned query sequence file. Fasta format",
              metavar="character"),
  make_option(c("-m", "--mummerprefix"), type="character", default="ref_query",
              help="Output name for MUMMER alignment files. [default= %default]",
              metavar="character"),
  make_option(c("-p", "--plot"), type="logical", default=TRUE,
              help="Plot corrected alignment. [default= %default]",
              metavar="logical"),
  make_option(c("-k", "--keep"), type="logical", default=TRUE,
              help="Keep alignments. [default= %default]",
              metavar="logical")
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

# Import packages ---------------------------------------------------
if (system('hostname',intern=T) == "taiji.embl.de") {
    .libPaths("~/R/x86_64-redhat-linux-gnu-library/3.2/")
}

library(ggplot2)
library(tidyr)
library(Biostrings)
library(optparse)


# run mummer
mummer_command = paste(paste("nucmer --prefix=", opt$mummerprefix, sep=""), "--coords", opt$ref, opt$query, sep = " ")
system(mummer_command, wait = TRUE)

# returns string w/o leading or trailing whitespace
trim <- function (x) {
  x = gsub("^\\s+|\\s+$", "", x)
  x = gsub("\\s+", "|", x)
  return(x)
}

# read in mummer table
mtable = read.delim(paste(opt$mummerprefix, ".coords", sep = ""),
  skip = 5, sep = "|", header = F, stringsAsFactors = F)
mtable = sapply(mtable, trim)
mtable = data.frame(mtable, stringsAsFactors = F)
mtable = mtable %>% separate(.,V1, c("S1","E1"), sep = "\\|") %>%
  separate(V2, c("S2","E2"), sep = "\\|") %>%
  separate(V3, c("LEN1","LEN2"), sep = "\\|") %>%
  separate(V5, c("REF","QUERY"), sep = "\\|")
# order table by ref start site
#mtable = mtable[order(as.numeric(mtable[,"S1"]), decreasing=F),]
mtable = mtable[order(as.numeric(mtable[,"LEN1"]), decreasing=T),]
mtable$S1 = as.numeric(mtable$S1)
mtable$E1 = as.numeric(mtable$E1)
mtable$S2 = as.numeric(mtable$S2)
mtable$E2 = as.numeric(mtable$E2)


# move query origin to match
qfa = readDNAStringSet(opt$query,"fasta")
qname = mtable[1,"QUERY"]
rname = mtable[1,"REF"]
qfa = qfa[[qname]]
#cat(as.numeric(mtable[1,"S2"]),"\n")
#cat(as.numeric(mtable[1,"E2"]),"\n")
if (mtable[1,"S2"] > mtable[1,"E2"]) {
  # match is reverse
  #cat("\n\n\nreverse_match\n\n\n\n\n")
  if (abs(mtable[1,"S1"]-mtable[1,"E2"]) <= 100) {
    # start is close - do nothing but reverse
    #cat("not moving origin\n")
    qfa_alt = reverseComplement(qfa)
  } else if (abs(mtable[1,"E2"] - 100) <= 100) {
    #cat("not moving origin\n")
    qfa_alt = reverseComplement(qfa)
    } else {
    #cat("moving origin\n")
    qfa_alt = reverseComplement(qfa)
    qfa_alt = DNAString(paste0(qfa[mtable[1,"E2"]:length(qfa)], qfa[1:(mtable[1,"E2"]-1)]))
  }
} else {
  #cat("\n\n\nforward_match\n\n\n\n\n")
  if (abs(mtable[1,"S1"]-mtable[1,"S2"]) <= 100) {
    # start is close - do nothing
    qfa_alt = qfa
  } else if (abs(mtable[1,"S2"] - 100) <= 100) {
    #cat("not moving origin\n")
    qfa_alt = qfa
  } else {
    qfa_alt = DNAString(paste0(qfa[mtable[1,"S2"]:length(qfa)], qfa[1:(mtable[1,"S2"]-1)]))
  }
}
if (is.null(opt$out)) {
  opt$out = paste(strsplit(opt$query, split = "\\.fasta")[[1]],
    "_",qname,"_",rname,"originAlign",".fasta",sep="")
}
write.fasta(qfa_alt, names = paste(qname, "_origin_fix", sep = ""),
  file.out = opt$out)

if (opt$plot) {
  # realign and plot
  ref = paste("REF=", opt$ref, sep="")
  query = paste("QUERY=", opt$out, sep="")
  prefix = paste("PREFIX=", opt$mummerprefix, "_originFix",sep="")
  mummer_command2 = paste(ref, query, prefix,
  "/usr/local/bin/MUMmer3.23/nucmer --maxmatch -c 100 -p $PREFIX $REF $QUERY",
  "/usr/local/bin/MUMmer3.23/mummerplot --postscript -p $PREFIX ${PREFIX}.delta -R $REF -Q $QUERY",
  "/usr/local/bin/MUMmer3.23/mummerplot --png -p $PREFIX ${PREFIX}.delta -R $REF -Q $QUERY",
  "ps2pdf ${PREFIX}.ps ${PREFIX}.pdf",
  sep=" && ")
  system(mummer_command2, wait = TRUE)
}

if (!opt$keep) {
  # remove everything expect plots
  rm_command = paste(
    "rm *.coords",
    "rm *.delta",
    sep = " && "
    )
  system(rm_command, wait = TRUE)
}
