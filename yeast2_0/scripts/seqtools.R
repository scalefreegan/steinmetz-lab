#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# Tools for tRNA-seq analysis
#
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2016"
.credits = ["Aaron Brooks"]
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"

# Import packages ---------------------------------------------------
library(ggplot2)
library(dplyr)
library(reshape2)
library(rtracklayer)
library(Rsamtools)
library(Biostrings)
library(seqinr)
library(robustbase)

remapRead = function(x) {
    #' Remap a read, trying to find the tRNA genomic
    #' location to which it maps best
    #'
    #' @param x A combined data.frame from GtRNA DB and bowtie2 alignment filtered for read
    #' @param returnBest Return the best match
    #' @return Alignment score
    #print(x)
    tRNA_strand = x$Strand
    read_strand = x$strand
    seq = as.character(x$seq)
    quality = as.character(x$qual)
    upstream = x$Upstream[[1]]
    downstream = x$Downstream[[1]]
    # do tRNA strand and mapped read strand agree?
    if (tRNA_strand != read_strand) {
        cat("tRNA and read strand do not agree")
    }
    #print(read_strand)
    if ((tRNA_strand == "+") && (read_strand == "+")) {
        # easy case
        # is read start closer to upstream or downstream region
        if (x$pos > x$Start) {
            d_name = "End"
        } else {
            d_name = "Start"
        }
        #print(d_name)
        if (d_name == "End") {
            # read should be aligned to downstream region
            # take only end region of read for alignment
            read_subseq = subseq(seq, abs(x$pos - x$End) + 2,nchar(seq))
            read_quality = subseq(quality, abs(x$pos - x$End) + 2,nchar(quality))
            pattern = DNAString(downstream)
            #print(seq)
            #print(DNAString(read_subseq))
            #print(pattern)

        } else if (d_name == "Start") {
            # read should be aligned to upstream region
            # take only start region of read for alignment
            read_subseq = subseq(seq, 1, abs(x$Start - x$pos))
            read_quality = subseq(quality, 1, abs(x$Start - x$pos))
            pattern = DNAString(upstream)
            #print(seq)
            #print(DNAString(read_subseq))
            #print(pattern)
        }
    } else if ((tRNA_strand == "-") && (read_strand == "-")) {
        # is read start closer to upstream or downstream region
        # note orientation is reversed due to limitations in reading bam file with Ramtools
        if (x$pos < x$Start) {
            d_name = "End"
        } else {
            d_name = "Start"
        }
        #print(d_name)
        if (d_name == "End") {
            # read should be aligned to downstream region
            # take only end region of read for alignment
            read_subseq = as.character(reverseComplement(DNAString(subseq(seq, 1, abs(x$pos - x$Start)))))
            read_quality = reverse(subseq(quality, 1, abs(x$pos - x$Start)))
            pattern = DNAString(downstream)
            #print(seq)
            #print(DNAString(read_subseq))
            #print(pattern)

        } else if (d_name == "Start") {
            # read should be aligned to upstream region
            # take only start region of read for alignment
            read_subseq = as.character(reverseComplement(DNAString(subseq(seq, abs(x$pos - x$End) + 2, nchar(seq)))))
            read_quality = reverse(subseq(quality, abs(x$pos - x$End) + 2, nchar(seq)))
            pattern = DNAString(upstream)
            #print(seq)
            #print(DNAString(read_subseq))
            #print(pattern)
        }
    }
    o = pairwiseAlignment(pattern = pattern, subject = DNAString(read_subseq),
                          subjectQuality = PhredQuality(read_quality),
                          gapOpening = 0, gapExtension = -5, type = "overlap",
                          scoreOnly = T)
    return(o)
}
