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
library(dplyr)
library(ggplot2)
library(seqinr)
library(Biostrings)
library(stringr)
library(parallel)
options("mc.cores"= 20)

tryload1 = try({
  f4 = "/g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scramble/SynIXR_segmentTable.rda"
  if (!file.exists(f4)) {
      segmentTable = str_split(as.character(JS94_alt_synIXR_corrected), loxPseq)[[1]]
      segmentTable = data.frame(number = seq(1:length(segmentTable)), seq = segmentTable)
      save(segmentTable, file = f4)
  } else {
      load(f4)
  }}
)
if (class(tryload1)=="try-error") {
  warning("Could not load segment table. You will need to make one manually")
}

tryload2 = try({
  scrambleTable = read.delim("/g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/synIXR/scramble_wpacbio.txt", sep="\t")
})
if (class(tryload2)=="try-error") {
  warning("Could not load scramble table. You will need to load it manually")
}

seg2seq = function(segmentOrder = c(1,2,-3), segmentTable, file = NA, sname = "") {
    #' Assemble sequence from scramble segment order. Optionally save to file as fasta
    #'
    #' @param segmentOrder Vector. Segment order. Negative values indicate inversion
    #' @param segmentTable Data.frame. Two columns: "number" segment integer, "seq", character
    #' @param file, Character. File path.
    #' @param sname, Character. Name of sequence for fasta header
    #' @return sequence
    #print(x)
    loxPseq = "ATAACTTCGTATAATGTACATTATACGAAGTTAT"
    loxPseq_rev = reverseComplement(DNAString(loxPseq))
    #loxPseq_rev = reverse(loxPseq)
    # assume loxPseq site after last base
    fseq = ""
    for (i in segmentOrder) {
        thisseq = as.character(filter(segmentTable, number == abs(i))$seq)
        if (i < 0) {
            if (i == -1) {
                # no loxP site afterwards
                #fseq = paste0(fseq, reverse(thisseq))
                fseq = paste0(fseq, reverseComplement(DNAString(thisseq)))
            } else {
                #inversion
                #fseq = paste0(fseq, reverse(thisseq), loxPseq_rev)
                fseq = paste0(fseq, reverseComplement(DNAString(thisseq)), loxPseq_rev)
            }

        } else {
            if (i == 44) {
                # no loxP site afterwards
                fseq = paste0(fseq, thisseq)
            } else {
                # normal
                fseq = paste0(fseq, thisseq, loxPseq)
            }
        }
    }
    if (!is.na(file)){
        write.fasta(fseq, names = sname, file.out = file)
    }
    fseq = DNAString(fseq)
    return(fseq)
}

seq2seg = function(seq, segmentTable) {
    #' Assemble segment order from scramble sequence.
    #'
    #' @param seq Character. scramble sequence
    #' @param segmentTable Data.frame. Two columns: "number" segment integer, "seq", character
    #' @return segment order
    seq = as.character(seq)
    # find loxP sites
    # split at loxP sites
    loxPseq = "ATAACTTCGTATAATGTACATTATACGAAGTTAT"
    lmatches = matchPattern(loxPseq, seq, max.mismatch = 3, with.indels = T)
    message(paste0("This sequence contains ", length(lmatches), " loxP sites"))
    thissegments = do.call(rbind, mclapply(seq(1,length(lmatches)+1),function(i){
        if (i == 1) {
            # first segment
            cat("first\n")
            thisseq = substr(seq, 1, start(lmatches[1])-1)
        } else if (i > length(lmatches)){
            # last segment
            cat("last\n")
            thisseq = substr(seq, end(lmatches[i-1])+1, nchar(seq))
        } else {
            thisseq = substr(seq, end(lmatches[i-1])+1, start(lmatches[i])-1)
        }
    align_f = pairwiseAlignment(DNAStringSet(segmentTable$seq), DNAString(thisseq), scoreOnly = T)
    align_r = pairwiseAlignment(reverseComplement(DNAStringSet(segmentTable$seq)),
                                DNAString(thisseq), scoreOnly = T)
    bestscores = c(max(align_f), max(align_r))
    fORr = which(bestscores == max(bestscores))
    if (fORr == 1) {
        # best match is forward
        segmentNumber = as.numeric(segmentTable$number[which(align_f == max(align_f))])
    } else {
        # best match is reverse
        segmentNumber = -as.numeric(segmentTable$number[which(align_r == max(align_r))])
    }
    return(segmentNumber)
    }))
return(thissegments)
}
