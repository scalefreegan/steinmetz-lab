#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# Tools for tRNA-seq analysis
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
library(ggplot2)
library(dplyr)
library(reshape2)
library(rtracklayer)
library(Rsamtools)
library(Biostrings)
library(seqinr)
library(robustbase)

makeHappyBam = function(bam, tRNA_annotations) {
    # alter bam data.frame from Rsamtools so that it can be combined with tRNA infos
    happy_bam = do.call(rbind, mclapply(seq(1,length(bam)), function(i){
        x = bam[[i]]
        x$seq = as.character(x$seq)
        x$qual = as.character(x$qual)
        x$Chr = rep(as.character(tRNA_annotations[i, "Chr"]), length(x$seq))
        x$Start = rep(as.numeric(tRNA_annotations[i, "Start"]), length(x$seq))
        x$End = rep(as.numeric(tRNA_annotations[i, "End"]), length(x$seq))
        x$AS = as.integer(x$tag$AS)
        if (is.null(x$XS)) {
          x$XS = rep(NA, length(x$seq))
        } else {
          x$XS = as.integer(x$tag$XS)
        }
        x$NM = as.integer(x$tag$NM)
        x$tag = NULL
        o = data.frame(x)
        return(o)
    }))
    return(happy_bam)
}

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
    attempt = try({
            if ((tRNA_strand == "+") && (read_strand == "+")) {
                # easy case
                # is read start closer to upstream or downstream region
                if (x$pos > x$Start) {
                    d_name = "End"
                } else {
                    d_name = "Start"
                }
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
        }, silent = T)
    if (class(attempt) == "try-error") {
        return(0)
    } else {
        o = pairwiseAlignment(pattern = pattern, subject = DNAString(read_subseq),
                          subjectQuality = PhredQuality(read_quality),
                          gapOpening = 0, gapExtension = -5, type = "overlap",
                          scoreOnly = T)
    }
    return(o)
}

filterReads = function(x) {
  #' Filters reads to ensure that at least one base from the read overlaps
  #' leader or trailer sequence
  #'
  #' @param x data.frame. Merged reads and tRNA info for anticodon
  #' @return data.frame with purely coding reads removed
  x$keep = apply(cbind(abs(x$Start - x$pos), abs(x$End - x$pos)), 1, min) < (nchar(as.character(x$seq)) - 1)
  x = filter(x, keep == T)
  x$keep = NULL
  return(x)
}

scoreAllReads = function(x, filter = FALSE) {
    #' Score all reads for anticodon set
    #'
    #' @param x data.frame. Merged reads and tRNA info for anticodon
    #' @return vector of scores
    if (filter) {
      x = filterReads(x)
    }
    reads = unique(levels(x$qname)[x$qname])
    o = mclapply(reads, function(read) {
        #print(read)
        thisreadset = filter(x, qname == read)
        #print(thisreadset)
        oo = lapply(seq(1:dim(thisreadset)[1]), function(i) {
            xx = thisreadset[i,]
            ooo = as.numeric(remapRead(xx))
            names(ooo) = paste(read, xx$Chr, xx$Start, sep=":")
            return(unlist(ooo))
        })
        if (sum(oo==0) == length(oo)) { # don't return if nothing mapped
            return(NULL)
        } else {
            return(oo)
        }
    })
    o = data.frame(score = unlist(o))
    r_names = strsplit(rownames(o), split = ":")[[1]]
    o$read = sapply(rownames(o), function(i) strsplit(i, split = ":")[[1]][1])
    o$chr = sapply(rownames(o), function(i) strsplit(i, split = ":")[[1]][2])
    o$tRNA.start = sapply(rownames(o), function(i) strsplit(i, split = ":")[[1]][3])
    rownames(o) = NULL
    return(o)
}

assignReads = function(mapped_reads, tRNA_annotations, anticodon = "AGC", useBam = FALSE) {
  #' Tries to assign reads that map to multiple tRNA locations to their best
  #' location using pariwise alignment of leader/trailer sequences
  #'
  #' @param mapped_reads data.frame. Mapped reads read in from bam file. Requires
  #' @param tRNA_annotations data.frame. tRNA_annotations from GtRNAdb
  #' @param anticodon. 3-letter anticodon specification
  #' @param useBam. Boolean. Use alignment score directly from bam file rather than local alignment
  #' @return vector of scores

  anticodon_tRNAs = filter(tRNA_annotations, Anticodon %in% anticodon)
  tRNA_read = merge(anticodon_tRNAs, mapped_reads, by = c("Chr","Start","End"))
  tRNA_read = filterReads(tRNA_read)
  mapping_scores = scoreAllReads(tRNA_read)
}

detectOutliers = function(x, a = -4, b = 3, c = 1.5, index = T) {
    #' Detect outliers in data with an extension of standard boxplot rule to
    #' allow for assymmetric data distributions
    #' Based on: http://www.r-bloggers.com/finding-outliers-in-numerical-data/
    #' for positive values of the medcouple MC, the adjusted boxplot rule’s nominal data range is:
    #' [Q1 – c * exp(a * MC) * IQD, Q3 + c * exp(b * MC) * IQD ]
    #' while for negative medcouple values, the nominal data range is:
    #' [Q1 – c * exp(-b * MC) * IQD, Q3 + c * exp(-a * MC) * IQD ]
    #'
    #' @param x Vector of numerical values
    #' @param a
    #' @param b
    #' @param c
    #' @param index Boolean. Whether to return the index (T) of value (F) of outliers
    #' @return Value or index of outliers

    # Compute the ‘medcouple’, a robust concept and estimator of skewness. The medcouple is defined
    # as a scaled median difference of the left and right half of distribution, and hence not based on the
    # third moment as the classical skewness.
    x = as.numeric(x)
    MC = mc(x, na.rm = T)
    Q1 = quantile(x, probs = 0.25)
    Q3 = quantile(x, probs = 0.75)
    IQD = Q3 - Q1
    if (MC >= 0) {
        o = c(Q1 - c * exp(a * MC) * IQD, Q3 + c * exp(b * MC) * IQD)

    } else {
        o = c(Q1 - c * exp(-b * MC) * IQD, Q3 + c * exp(-a * MC) * IQD)
    }
    return(o)
}
