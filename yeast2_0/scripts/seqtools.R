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
try({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(rtracklayer)
  library(Rsamtools)
  library(Biostrings)
  library(seqinr)
  library(robustbase)
})

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
                #print(d_name)
                if (d_name == "End") {
                    # read should be aligned to downstream region
                    # take only end region of read for alignment
                    read_subseq = DNAString(subseq(seq, abs(x$pos - x$End) + 2,nchar(seq)))
                    read_quality = subseq(quality, abs(x$pos - x$End) + 2,nchar(quality))
                    pattern = DNAString(downstream)
                    pattern = pattern[1:length(read_subseq)]
                    #print(x$Chr)
                    #print(seq)
                    #print(DNAString(read_subseq))
                    #print(downstream)
                    #print(pattern)
                } else if (d_name == "Start") {
                    # read should be aligned to upstream region
                    # take only start region of read for alignment
                    read_subseq = DNAString(subseq(seq, 1, abs(x$Start - x$pos)))
                    read_quality = subseq(quality, 1, abs(x$Start - x$pos))
                    pattern = DNAString(upstream)
                    pattern = pattern[((length(pattern)-length(read_subseq))+1):length(pattern)]
                    #print(x$Chr)
                    #print(seq)
                    #print(DNAString(read_subseq))
                    #print(upstream)
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
                        read_subseq = DNAString(as.character(reverseComplement(DNAString(subseq(seq, 1, abs(x$pos - x$Start))))))
                        read_quality = reverse(subseq(quality, 1, abs(x$pos - x$Start)))
                        pattern = DNAString(downstream)
                        pattern = pattern[1:length(read_subseq)]
                        #print(x$Chr)
                        #print(seq)
                        #print(DNAString(read_subseq))
                        #print(downstream)
                        #print(pattern)

                    } else if (d_name == "Start") {
                        # read should be aligned to upstream region
                        # take only start region of read for alignment
                        read_subseq = DNAString(as.character(reverseComplement(DNAString(subseq(seq, abs(x$pos - x$End) + 2, nchar(seq))))))
                        read_quality = reverse(subseq(quality, abs(x$pos - x$End) + 2, nchar(seq)))
                        pattern = DNAString(upstream)
                        pattern = pattern[((length(pattern)-length(read_subseq))+1):length(pattern)]
                        #print(x$Chr)
                        #print(seq)
                        #print(DNAString(read_subseq))
                        #print(upstream)
                        #print(pattern)
                    }
            }
        }, silent = T)
    if (class(attempt) == "try-error") {
        return(0)
    } else {
        mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
        o = pairwiseAlignment(pattern = pattern, subject = read_subseq,
                          subjectQuality = PhredQuality(read_quality),
                          #substitutionMatrix = mat,
                          gapOpening = 0, gapExtension = -5, type = "local",
                          scoreOnly = T)
    }
    return(o)
}

filterMature = function(x, flip = FALSE, clean = TRUE) {
  #' Filters reads to ensure that at least one base from the read overlaps
  #' leader or trailer sequence. Categorize those that do not as likely
  #' mature. Additionaly, categorizes those that have a trailer sequence
  #' C, CC, CCA, or - on the other strand - G, G, GGA only as likely
  #' mature sequences
  #'
  #' @param x data.frame. Merged reads and tRNA info for anticodon
  #' @return data.frame with purely coding reads removed

  substrRight <- function(z, n) {
    substr(z, nchar(z)-n+1, nchar(z))
  }

  x$mature = sapply(seq(1,dim(x)[1]), function(i){
    i = x[i,]
    if (i$pos > i$Start) {
        d_name = "End"
    } else {
        d_name = "Start"
    }
    #print(d_name)
    if (i$strand == "+") {
        if (d_name == "Start") {
            i$toend = i$Start - i$pos
            i$mature = i$toend <= 0
        } else {
            i$toend = i$pos + nchar(as.character(i$seq)) - 1 - i$End
            # check for CCA
            i$mature = i$toend <= 0
            if (i$toend <= 3 & i$mature != TRUE) {
                lastchar = substrRight(as.character(i$seq), i$toend)
                if (lastchar %in% c("C","CC","CCA")) {
                    i$mature = TRUE
                } else {
                    i$mature = FALSE
                }
            }
        }
    } else if (i$strand == "-"){
        if (d_name == "Start") {
            i$toend = i$Start - i$pos
            i$mature = i$toend <= 0
            if (i$toend <= 3 & i$mature != TRUE) {
                # check for TGG
                lastchar = substr(as.character(i$seq), 1, i$toend)
                if (lastchar %in% c("T","TG","TGG")) {
                    i$mature = TRUE
                } else {
                    i$mature = FALSE
                }
            }
        } else {
            i$toend = i$pos + nchar(as.character(i$seq)) - 1 - i$End
            i$mature = i$toend <= 0
        }
    }
    #print(i)
    return(i$mature)
  })



  if (flip) {
    x = filter(x, mature == T)
  } else {
    x = filter(x, mature == F)
  }
  if (clean) {
    x$mature = NULL
  }
  return(x)
}

scoreAllReads = function(x, filter = FALSE, ...) {
    #' Score all reads for anticodon set
    #'
    #' @param x data.frame. Merged reads and tRNA info for anticodon
    #' @return vector of scores
    if (filter) {
      x = filterMature(x, ...)
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

assignReads = function(mapped_reads, tRNA_annotations, anticodon = "AGC", useBam = FALSE,
  returnCounts = FALSE, ...) {
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
    tRNA_read = filterMature(tRNA_read, ...)
    if (useBam) {
      mapping_scores = tRNA_read %>% group_by(qname) %>% do({
        maxs = max(.$AS)
        o = filter(., AS >= maxs)
        return(o)
      })
      mapped = data.frame(score = mapping_scores$AS, read = mapping_scores$qname,
        chr = unlist(mapping_scores$Chr), tRNA.start = unlist(mapping_scores$Start))
    } else {
      mapping_scores = scoreAllReads(tRNA_read)
      mapped = mapping_scores %>% group_by(read) %>% do({
        maxs = max(.$score)
        o = filter(., score >= maxs)
        return(o)
      })
    }
  uniquemap = round(sum(table(mapped$read)==1)/length(unique(mapped$read)) * 100, 2)
  cat(paste0(uniquemap, "% of reads mapped uniquely\n"))
  if (returnCounts) {
    mapped = mapped %>% group_by(chr, tRNA.start) %>% do({data.frame(counts = dim(.)[1])})
  }
  return(mapped)
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

processDAZZstats = function(dbstats) {
  # number of reads
  nreads = trimws(as.character(dbstats[1,1]))
  nreads = as.numeric(paste(strsplit(strsplit(nreads,split = " ")[[1]][1], split=",")[[1]], collapse=""))

  # average length
  avglength = strsplit(trimws(as.character(dbstats[3,1])), split = " ")[[1]][[1]]
  avglength = as.numeric(paste(strsplit(avglength, split=",")[[1]], collapse = ""))

  # sd length
  sdlength = strsplit(trimws(as.character(dbstats[4,1])), split = " ")[[1]][[1]]
  sdlength = as.numeric(paste(strsplit(sdlength, split=",")[[1]], collapse = ""))

  # base composition
  basecomp = trimws(strsplit(trimws(as.character(dbstats[5,1])), split = ":")[[1]][[2]])
  basecomp = strsplit(basecomp, split = " ")[[1]]
  basecomp = sapply(basecomp, function(i){
      tor = strsplit(i, split = "\\(")[[1]]
      tor[2] = sub(pattern = "\\)", replacement = "",  tor[2])
      tor = rev(tor)
      return(tor)
  })
  colnames(basecomp) = basecomp[1,]
  tmpbasecomp = as.numeric(basecomp[2,])
  names(tmpbasecomp) = colnames(basecomp)
  basecomp = tmpbasecomp

  # base distribution, histogram
  basecomp_header = trimws(strsplit(trimws(as.character(dbstats[7,1])), split = ":")[[1]])
  basecomp_header[2] = gsub("%","p",basecomp_header[2])
  basecomp_header[2] = strsplit(basecomp_header[2], split = " ")
  basecomp_header[[2]] = basecomp_header[[2]][basecomp_header[[2]] != ""]
  basecomp_header[[2]] = sapply(list(c(1,1),c(2,3),c(4,5),c(6,6)), function(i){
      #print(i)
      paste(basecomp_header[[2]][i[1]:i[2]], collapse="")
  })
  basedist = data.frame(dbstats[8:nrow(dbstats),1], stringsAsFactors = F)
  colnames(basedist) = "start"
  basedist = basedist %>% rowwise() %>%
      do({
          vals = gsub("\\s+", " ", str_trim(as.character(.$start)))
          vals = sub(": ", ":", vals)
          #print(levels(basedist$start)[.])
          #print(as.character(.$start))
          data.frame(start = vals, stringsAsFactors = F)
      })
  basedist = basedist %>% separate(col = start, into = c("length","rest"),sep = ":") %>%
      separate(col = rest, into = c("a","b", "c", "d"),sep = " ")
  colnames(basedist) = unlist(basecomp_header)
  basedist$Bin = as.numeric(gsub(",","",basedist$Bin))
  basedist$Count = as.numeric(gsub(",","",basedist$Count))
  basedist = apply(basedist,2,as.numeric)
  o = list()
  o$nreads = nreads
  o$avglength = avglength
  o$sdlength = sdlength
  o$basecomp = basecomp
  o$basedist = basedist
  return(o)
}
