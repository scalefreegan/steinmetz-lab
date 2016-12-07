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
      segmentTable$essential = F
      segmentTable$essential[c(7,9,10,12,20)] = T
      segmentTable = select(segmentTable,number,essential,seq)
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

# clean up
rm("tryload1","tryload2","f4")

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

seq2seg = function(seq, segmentTable, segThreshold = 0.33) {
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
            # cat("first\n")
            thisseq = substr(seq, 1, start(lmatches[1])-1)
        } else if (i > length(lmatches)){
            # last segment
            # cat("last\n")
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
            bmatch = which(align_f == max(align_f))
        } else {
            # best match is reverse
            bmatch = which(align_r == max(align_r))
        }
        if (nchar(thisseq) <= (segThreshold * nchar(as.character(segmentTable$seq[bmatch])))) {
          segmentNumber = NULL
        } else {
          if (fORr == 1) {
            segmentNumber = as.numeric(segmentTable$number[bmatch])
          } else {
            segmentNumber = -as.numeric(segmentTable$number[bmatch])
          }
        }
        return(segmentNumber)
    }))
return(as.vector(thissegments))
}

drawSegments = function(segments = c(1,2,3), thissegmentTable = segmentTable, plot = T, lwd = .5,
  width = 5, height = 2, label = TRUE, dev.control = FALSE, label.size = 5, label.y = 0) {
    segmentDF = data.frame(number = abs(segments), raw = segments, order = seq(1,length(segments)))
    subTable = merge(segmentDF,thissegmentTable, by = "number")
    plotDF = subTable %>% group_by(order) %>% do({
        thisdf = .[1,]
        i = as.numeric(thisdf$order)
        if (thisdf$raw > 0) {
            o = data.frame(number=rep(thisdf$number,3),raw=rep(thisdf$raw,3),
                   polygon.x = c(i-1,i-1,i), polygon.y = c(0,2,1),
                   essential=rep(thisdf$essential,3))
        } else {
            o = data.frame(number=rep(thisdf$number,3),raw=rep(thisdf$raw,3),
                   polygon.x = c(i,i,i-1), polygon.y = c(0,2,1),
                   essential=rep(thisdf$essential,3))
        }
        return(o)
    })
    p <- ggplot(data = plotDF)
    xpos = plotDF %>% group_by(order) %>% summarise(x = mean(polygon.x))
    ypos = plotDF %>% group_by(order) %>% summarise(y = label.y)
    thislabel = plotDF %>% group_by(order) %>% summarise(label = unique(number))
    p <- p + geom_polygon(aes(x=polygon.x,y=polygon.y,group=order,color=essential,fill=number),size=lwd) +
      xlab("") + ylab("") +
      theme(axis.text = element_text(size = 0),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white"),
            legend.position = "none",
            axis.ticks = element_blank()) +
      scale_fill_gradientn(colours = rev(rainbow(dim(thissegmentTable)[1],start=6/6,end=5/6)),limits=c(0,dim(thissegmentTable)[1])) +
      scale_colour_manual(values=c("white","red"))
    if (label) {
      p <- p + annotate("text", x = xpos$x, y = ypos$y, label = thislabel$label, size = label.size)
    }
    if (dev.control) {
      dev.new(width=width, height=height)
    }
    p
}

loc2seg = function(start, end = NA, segment_order = NA, sequence = NA,
  thissegmentTable = segmentTable, concise = F, force_forward = F, ...) {
  if (is.na(segment_order) && is.na(sequence)) {
    cat("ERROR: Must supply either segment_order or sequence")
    return(NULL)
  }
  if (!is.na(sequence)) {
    # convert sequence to segment order
    segment_order = seq2seg(sequence, thissegmentTable, ...)
  }
  df = data.frame(number = segment_order, order = seq(1,length(segment_order)))
  df = merge(df, thissegmentTable, by = "number", keep.x = T) %>% arrange(order)
  # make sure it's a character. should make this a global fix
  df$seq = as.character(df$seq)
  if (!is.na(end)) {
  # always work with coordinates on Watson (+) strand
    if (!force_forward) {
      if (start > end) {
        internal.start = end
        internal.end = start
      } else {
        internal.start = start
        internal.end = end
      }
    }
  } else {
    internal.start = start
    internal.end = start
  }
  # set up vectors for addtions
  order.v = df$order
  start.v = rep(0,length.out = length(order.v))
  end.v = rep(0,length.out = length(order.v))
  width.v = rep(0,length.out = length(order.v))
  for (i in 1:length(df$order)) {
      #print(i)
      order.v[i] = df$order[i]
      thisrow = df %>% filter(order==df$order[i])
      if (order.v[i] == 1) {
        start.v[i] = 1
        width.v[i] = as.numeric(nchar(thisrow$seq))
        end.v[i] = width.v[i]

      } else {
        start.v[i] = end.v[i-1] + 1
        width.v[i] = as.numeric(nchar(thisrow$seq))
        end.v[i] = start.v[i] + (width.v[i] - 1)
      }
  }
  df = data.frame(df,start = start.v,end = end.v,width = width.v) %>%
    select(order, number, essential, start, end, width, seq)
  tor.s = df %>% filter(internal.start >= start, internal.start <=  end)
  tor.e = df %>% filter(internal.end >= start, internal.end <= end)
  if (!force_forward) {
    s.e = sort(c(min(tor.s$order), max(tor.e$order)))
    tor.o = seq(s.e[1], s.e[2])
    tor = df %>% filter(order %in% tor.o)
  } else {
    tor.s = arrange(tor.s, start)
    tor.s = data.frame(force.order = seq(1,length.out = dim(tor.s)[1]),tor.s)
    tor.e = arrange(tor.e, start)
    tor.e = data.frame(force.order = seq((max(tor.s$force.order)+1),length.out = dim(tor.e)[1]),tor.e)
    tor = bind_rows(tor.s, tor.e) %>% arrange(abs(start-internal.start))
  }
  if (concise) {
    if (force_forward) {
      return(paste(tor$number,collapse="-"))
    }
    else if (start > end) {
      tor = arrange(tor,-order)
      return(paste(tor$number,collapse="-"))
    } else {
      tor = arrange(tor,order)
      return(paste(tor$number,collapse="-"))
    }
  } else {
    tor = data.frame(query.start = start, query.end = end, tor)
    return(tor)
  }
}
