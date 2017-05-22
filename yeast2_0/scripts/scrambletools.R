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
library(lsa)
library(kernlab)
options("mc.cores"= 20)

tryload1 = try({
  f4 = "/g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scramble/SynIXR_segmentTable.rda"
  f4_txt = "/g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scramble/SynIXR_segmentTable.txt"
  if (!file.exists(f4)) {
      segmentTable = str_split(as.character(JS94_alt_synIXR_corrected), loxPseq)[[1]]
      segmentTable = data.frame(number = seq(1:length(segmentTable)), seq = segmentTable)
      segmentTable$essential = F
      segmentTable$essential[c(7,9,10,12,20)] = T
      segmentTable = select(segmentTable,number,essential,seq)
      save(segmentTable, file = f4)
      write_delim(segmentTable, delim="/t", path = f4_txt, col_names = T)
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

seg2seq <- function(segmentOrder = c(1,2,-3), segmentTable, file = NA, sname = "") {
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

seq2seg <- function(seq, segmentTable, segThreshold = 0.33) {
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

drawSegments <- function(segments = c(1,2,3), thissegmentTable = segmentTable, plot = T, lwd = .5,
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
    #print(head(plotDF))
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

loc2seg <- function(start, end = NA, segment_order = NA, sequence = NA,
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

strain2seg <- function(strain = NA, scrambleTable = NA) {
  if (is.na(strain) || is.na(scrambleTable)) {
    cat("ERROR: Must supply either strain and scrambleTable")
    return(NULL)
  }
  thisstrain = scrambleTable %>% filter(Strain.ID %in% strain)
  tor = lapply( seq(1,dim(thisstrain)[1]), function(i){
      thisseg = as.character(thisstrain$solutions[[i]])
      thisseg = strsplit(thisseg,split="\n")[[1]]
      torr = lapply(thisseg, function(j){
          if (grepl(",",j)==F) {
            j = strsplit(j, split="-")[[1]]
            #print(j)
            j = seq(as.numeric(j[1]),as.numeric(j[2]))
          }
          else {
            j = as.numeric(strsplit(j,split=",")[[1]])
          }
          return(j)
      })
      return(torr)
  })
  names(tor) = thisstrain$Strain.ID
  return(tor)
}

readgff <- function(f) {
  gffnames=c('seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute')
  gff = read_delim(f,delim="\t", col_names = gffnames)
  gff = gff %>%
        mutate(ID = sub(";.*$","",sub("^ID=","",attribute)),parent = sub("-[0-9]$","",ID)) %>%
        select(seqname, source, feature, start, end, score, strand, ID, parent, attribute)
  return(gff)
}

readbed <- function(f) {
  bed_cnames = c(
    "sample",
    "seqname",
    "start",
    "end",
    "id",
    "score",
    "strand",
    "cigar"
    )
  bed = read_delim(f,delim = "\t",col_names = bed_cnames)
  return(bed)
}

getAllCopies <- function(segment, gff) {
  gff %>% filter(parent == segment)
}

getOverlaps <- function(bed, gffline) {
  gff_seqname = gffline$seqname
  gff_start = gffline$start
  gff_end = gffline$end
  thisbed_f = bed %>%
            filter(seqname %in% gff_seqname, end >= gff_start)
  thisbed_r = bed %>%
            filter(seqname %in% gff_seqname, start <= gff_end)
  out = dplyr::intersect(thisbed_f, thisbed_r)
  return(out)
}

getCounts <- function(bed, gffline, crop=T, multicore=F) {
  tor = bed %>% group_by(sample) %>% do({
    thismin = min(.$start)
    thismax =  max(.$end)
    out = rep(0,length(seq(1,(thismax-thismin+1))))
    names(out) = seq(thismin,thismax)
    #print(thismin)
    #print(thismax)
    if (multicore) {
      x = mclapply(seq(1,dim(.)[1]),function(i){
        #print(i)
        thisrow = .[i,]
        theseinds = seq(thisrow$start,thisrow$end)
        #out[theseinds] = out[theseinds] + 1
        })
    } else {
      x = lapply(seq(1,dim(.)[1]),function(i){
        #print(i)
        thisrow = .[i,]
        theseinds = seq(thisrow$start,thisrow$end)
        #out[theseinds] = out[theseinds] + 1
        })
    }
    thesecounts = table(unlist(x))
    if (gffline$strand == "-") {
      thesecounts = thesecounts[order(names(thesecounts))]
      names(thesecounts) = rev(names(thesecounts))
    }
    out[names(thesecounts)] = thesecounts[names(thesecounts)]
    out = data.frame(position = as.numeric(names(out)), count = out)
    out = arrange(out,position)
    out = out %>%
          mutate(relative_pos = position - gffline$start,
                  ID = paste(gffline$seqname,gffline$ID,sep="/"),
                  copy = ifelse(grepl("-",gffline$ID),gsub("^.*-","",gffline$ID),as.character(1)))
    if (crop) {
      # only return counts from this segment, ignore small adjacent region_start
      out = out %>% filter(relative_pos >= 0 & relative_pos <= (gffline$end-gffline$start))
    }
    return(out)
  })
  return(tor)
}

getCounts_multi <- function(bed, gff_one_seg, gff) {
  thisseg = as.character(unique(gff_one_seg$parent))
  other_gffs = getAllCopies(thisseg, gff)
  o = lapply(seq(1,dim(other_gffs)[1]), function(x){
      #print(x)
      thisline = other_gffs[x,]
      thisovlps = getOverlaps(bed,thisline)
      if (dim(thisovlps)[1]>0) {
        thiscounts = getCounts(thisovlps,thisline,crop=T)
        #print(head(thiscounts))
        return(thiscounts)
      } else {
        return(data.frame())
      }
    })
  # get rid of empty data.frames
  o = o[sapply(o,function(x)dim(x)[1]>0)]
  if (length(o)==0) {
    return(data.frame())
  } else {
    o = do.call(rbind,o)
    o = o %>% group_by(sample, relative_pos) %>%
            summarise(count=sum(count)) %>%
            arrange(relative_pos) %>%
            mutate(segment = thisseg,
                   ID = gff_one_seg$ID,
                   seqname = unique(gff_one_seg$seqname))
    return(o)
  }
}

translateFromRelative <- function(ID, relative_pos,thisgff) {
  start = sapply(ID, function(i) {thisgff$start[which(thisgff$ID == i)]})
  df = data.frame(ID,relative_pos, start)
  df = df %>% mutate(abs_position = relative_pos + start)
  abs_position = df$abs_position
  return(abs_position)
}

getCounts_full <- function(targetSegment, bed, thisgff, gff) {
  target_gff = thisgff %>% filter(ID == targetSegment)
  #print(target_gff)
  # target_ovlps = getOverlaps(bed,target_gff)
  # target_counts = getCounts(target_ovlps, target_gff,multicore=F)
  target_counts = getCounts_multi(bed,target_gff,gff)
  if (dim(target_counts)[1]==0) {
    return(data.frame())
  }
  #print(head(target_counts))
  target_counts = target_counts %>% mutate(segment=target_gff$parent)
  other_gff = thisgff %>% filter(ID != targetSegment)
  other_gff = other_gff %>% filter(feature=="engineered_region")
  if (dim(other_gff)[1]>0) {
    other_counts = lapply(seq(1,dim(other_gff)[1]), function(i){
          tor = getCounts_multi(bed,other_gff[i,],gff)
          return(tor)
        })
        # get rid of empty data.frames
        other_counts = other_counts[sapply(other_counts,function(x)dim(x)[1]>0)]
        if (length(other_counts)==0) {
        return(data.frame())
    } else {
        other_counts = do.call(rbind,other_counts)
    }
    full_counts = rbind(target_counts,other_counts)
  } else {
    full_counts = target_counts
  }
  segOrder = thisgff %>% filter(feature=="engineered_region") %>% select(parent)
  segOrder = segOrder[[1]]
  thissegment = gsub("-.*$","",targetSegment)
  full_counts = full_counts %>%
                mutate(segment=factor(segment,levels=segOrder)) %>%
                arrange(segment,relative_pos)
  full_counts = full_counts %>% group_by(sample,segment) %>%
                 mutate(position = translateFromRelative(ID, relative_pos,thisgff))
  full_counts = full_counts %>%
                mutate(relative_pos = position - target_gff$start,
                       copy = ifelse(segment==thissegment & grepl("-",targetSegment),
                                    gsub("^.*-","",targetSegment),as.character(1)),
                        targetSegment=targetSegment) %>%
                select(sample,segment,seqname,position,relative_pos,count,copy,targetSegment,ID) %>% ungroup()
  return(full_counts)
}

getAdjacentSegs <- function(segment, gff, adjacent) {
  # make sure gff is sorted
  gff_segment = gff %>%
                filter(feature=="engineered_region") %>%
                arrange(seqname,start)
  theseinds = which(gff_segment$parent == segment)
  seg_gff = gff_segment[theseinds,]
  thisgff = lapply(theseinds, function(i){
    if (adjacent > 0) {
      o = gff_segment[seq(i-adjacent,i+adjacent),]
    } else {
      o = gff_segment[i,]
    }
    region_start = min(o$start)
    region_end = max(o$end)
    gff_f = gff %>%
      filter(feature!="engineered_region",
            seqname %in% as.character(unique(o$seqname)),
            end >= region_start)
    gff_r = gff %>% filter(feature!="engineered_region",
            seqname %in% as.character(unique(o$seqname)),
            start <= region_end)
    gff_gene = dplyr::intersect(gff_f, gff_r)
    o = rbind(o,gff_gene)
    return(o)
  })
  names(thisgff) = seg_gff$ID
  return(thisgff)
}

aggCount <- function(segment,bed,gff,seqnames=NA,adjacent=2) {
  thisgff = getAdjacentSegs(segment,gff,adjacent)
  if (length(thisgff)==0) {
    o = list()
    o[["gff"]] = list()
    o[["count"]] = data.frame()
  } else {
    thiscounts = do.call(rbind,lapply(names(thisgff), function(i){
      #print(i)
      getCounts_full(i,bed,thisgff[[i]],gff)
    }))
    thiscounts = thiscounts %>% mutate(query=segment)
    # change segment to numeric
    thiscounts$segment = as.numeric(levels(thiscounts$segment)[thiscounts$segment])
    seqname = sapply(thisgff, function(i)unique(i$seqname))
    names(thisgff) = paste(seqname,names(thisgff),sep=":")
    o = list()
    o[["gff"]] = thisgff
    o[["count"]] = thiscounts

  }
  return(o)
}

aggCount_multigff_multiseg <- function(segments,bed,gff_list,adjacent) {
  thisdf = mclapply(segments, function(thisseg){
      print(thisseg)
      thesecounts = lapply(seq(1,length(gff_list)),function(thisgff){
          x = aggCount(thisseg,bed,gff_list[[thisgff]],adjacent=adjacent)
          # this is specific to naming scheme of current
          #x$count = x$count %>% separate(sample,into=c("strain","rep"),sep = "x",remove = F)
        })
      names(thesecounts) = names(gff_list)
      o = list()
      o$count
      o$gff
      for (i in names(gff_list)) {
        o$count = rbind(o$count,thesecounts[[i]][["count"]])
        o$gff = c(o$gff,thesecounts[[i]][["gff"]])
      }
      return(o)
    })
  names(thisdf) = segments
  return(thisdf)
}

cosine <- function (x, y)
{
    if (is.matrix(x) && is.null(y)) {
        co = array(0, c(ncol(x), ncol(x)))
        f = colnames(x)
        dimnames(co) = list(f, f)
        for (i in 2:ncol(x)) {
            for (j in 1:(i - 1)) {
                co[i, j] = cosine(x[, i], x[, j])
            }
        }
        co = co + t(co)
        diag(co) = 1
        return(as.matrix(co))
    }
    else if (is.vector(x) && is.vector(y)) {
        return(crossprod(x, y)/sqrt(crossprod(x) * crossprod(y)))
    }
    else if (is.vector(x) && is.matrix(y)) {
        co = vector(mode = "numeric", length = ncol(y))
        names(co) = colnames(y)
        for (i in 1:ncol(y)) {
            co[i] = cosine(x, y[, i])
        }
        return(co)
    }
    else {
        stop("argument mismatch. Either one matrix or two vectors needed as input.")
    }
}
class(cosine) <- "kernel"

makeKernel <- function(thisdf,bed,idealKernel=F,long=F,log=F) {
  allsamples = unique(bed$sample)
  if (long) {
    mat = lapply(seq(1,length(thisdf)),function(i){
      #print(i)
      thisdf[[i]]$sample = factor(thisdf[[i]]$sample, levels=allsamples)
      thismat = acast(data = thisdf[[i]]$count,
                  formula = factor(sample,levels=allsamples) ~ relative_pos,
                  value.var = "count",
                  fun.aggregate = mean,
                  fill=0)
      if (log) {
        print("doing log transform")
        thismat = log10(thismat+1)
      }
      idealmat = matrix(data = 0,
        nrow = length(allsamples),
        ncol=dim(thismat)[2],
        dimnames = list(allsamples,colnames(thismat)))
      idealmat[rownames(thismat),colnames(thismat)] = thismat[rownames(thismat),colnames(thismat)]
      return(idealmat)
    })
    mat = do.call(cbind,mat)
    k = kernelMatrix(cosine,mat)
    rownames(k) = allsamples
    colnames(k) = allsamples
    o = list()
    o$mat = mat
    o$kernel = k
  } else {
    mat = acast(data = thisdf$count,
                formula = sample ~ relative_pos,
                value.var = "count",
                fun.aggregate = mean,
                fill=0)
    if (log) {
      #print("doing log transform")
      mat = log10(mat+1)
    }
    k = kernelMatrix(cosine,mat)
    rownames(k) = rownames(mat)
    colnames(k) = rownames(mat)
    o = list()
    if (idealKernel) {
      # return a kernel with all samples represented
      idealk = matrix(data = 1,
        nrow = length(allsamples),
        ncol=length(allsamples),
        dimnames = list(allsamples,allsamples))
      idealk[rownames(k),colnames(k)] = k[rownames(k),colnames(k)]
      newentries = setdiff(allsamples,rownames(k))
      toremove = expand.grid(newentries,rownames(k))
      #print(toremove)
      idealk[as.vector(toremove[,1]),as.vector(toremove[,2])] = 0
      idealk[as.vector(toremove[,2]),as.vector(toremove[,1])] = 0
      o$kernel = idealk
    } else {
      o$kernel = k
    }
    o$mat = mat
  }
  return(o)
}

combineKernels <- function(kernelList) {
  k = kernelList[[1]][["kernel"]]
  if (length(kernelList)>1) {
    for (i in 2:length(kernelList)) {
      k = k + kernelList[[i]][["kernel"]]
    }
  }
  k = as.kernelMatrix(k)
  return(k)
}

plotCounts <- function(thisdf, thisstart=NA,thisend=NA,
  segmentTable=segmentTable, rect_h = 0.1,
  rect_s = 0.15, shift=F, logy = F, useReps = T, plot_genes = F) {
  segmentTable = segmentTable %>%
                  mutate(color = rev(rainbow(dim(segmentTable)[1],start=6/6,end=5/6))[number])
  maxcount = max(thisdf$count$count)
  maxcopy = max(thisdf$count$copy)
  rectdf = thisdf$count %>% group_by(sample,targetSegment,segment,strain) %>% do({
    thisgroup = .
    thistargetSegment = as.character(unique(thisgroup$targetSegment))
    thissegment = as.character(unique(thisgroup$segment))
    thisseqname = as.character(unique(thisgroup$seqname))
    thisgff = thisdf$gff[[paste(thisseqname,thistargetSegment,sep=":")]]
    target_gff = thisgff %>% filter(ID==thistargetSegment)
    gff_entry = thisgff %>% filter(parent==thissegment)
    # convert to relative coordinates
    xmin = gff_entry$start - target_gff$start
    xmax = gff_entry$end - target_gff$start
    if (logy) {
      ymin = maxcount + (maxcount * (rect_s))
      yrat1 = ymin/maxcount
      ymin = 10^(log10(maxcount)*yrat1)
      ymax = ymin + (maxcount * rect_h)
      yrat2 = ymax/ymin
      ymax = 10^(log10(ymin)*yrat2)
    } else {
      ymin = maxcount + (maxcount * (rect_s))
      ymax = ymin + (maxcount * rect_h)
    }
    #print(thisgroup$segment)
    #color = segmentTable[which(segmentTable$number==unique(thisgroup$segment)),"color"]
    tor = data.frame(xmin,xmax,ymin,ymax)
    return(tor)
  })
  # if (shift) {
  #   # adjust y values based on copy number
  #   # only when plotting together
  #   rectdf = rectdf %>% group_by(sample) %>% do({
  #     thisgroup = .
  #     these_targetSegment = as.character(unique(thisgroup$targetSegment))
  #     if (length(these_targetSegment)>1) {
  #       o = thisgroup
  #       for (i in 2:length(these_targetSegment)) {
  #         previous = filter(o, targetSegment==these_targetSegment[i-1])
  #         previous_ymin = as.numeric(unique(previous$ymin))
  #         previous_ymax = as.numeric(unique(previous$ymax))
  #         if (logy) {
  #           new_ymin = previous_ymax + (previous_ymax *rect_s)
  #           yrat1 = new_ymin/previous_ymax
  #           new_ymin = 10^(log10(previous_ymax)*yrat1)
  #           new_ymax = new_ymin + (maxcount * rect_h)
  #           yrat2 = new_ymax/new_ymin
  #           new_ymax = 10^(log10(new_ymin)*yrat2)
  #         } else {
  #           new_ymin = previous_ymax + (previous_ymax *rect_s)
  #           new_ymax = new_ymin + (maxcount * rect_h)
  #         }
  #         o = o %>% mutate(ymin = ifelse(targetSegment==these_targetSegment[i],new_ymin,ymin),
  #                         ymax = ifelse(targetSegment==these_targetSegment[i],new_ymax,ymax))
  #       }
  #     } else {
  #       o = thisgroup
  #     }
  #     return(o)
  #     })
  # }

  polygon_dfs = thisdf$count %>% group_by(sample,targetSegment,segment,strain) %>% do({
    thisgroup = .
    thistargetSegment = as.character(unique(thisgroup$targetSegment))
    thissegment = as.character(unique(thisgroup$segment))
    thisseqname = as.character(unique(thisgroup$seqname))
    thisgff = thisdf$gff[[paste(thisseqname,thistargetSegment,sep=":")]]
    target_gff = thisgff %>% filter(ID==thistargetSegment)
    gff_entry = thisgff %>% filter(parent==thissegment)
    # make gene polygons
    absStart = as.numeric(unique(filter(thisgff,ID==thistargetSegment) %>% select(start)))
    # cheating here
    ymin = min(rectdf$ymin)
    ymax = min(rectdf$ymax)
    ymid = mean(c(ymin,ymax))
    #print(ymid)
    features2draw = thisgff %>% filter(feature != "engineered_region") %>%
                    mutate(relative_start=start-absStart,
                                relative_end=end-absStart)
    if (dim(features2draw)[1]>0) {
      o = do.call(rbind,lapply(seq(1,dim(features2draw)[1]),function(i){
          #print(i)
          #print(as.numeric(features2draw[i,"relative_start"]))
        if (features2draw[i,"strand"]=="+") {
          polygon_x = c(as.numeric(features2draw[i,"relative_start"]),
                        as.numeric(features2draw[i,"relative_start"]),
                        as.numeric(features2draw[i,"relative_end"]))
          polygon_y = c(ymin,ymax,ymid)
        } else {
          polygon_x = c(as.numeric(features2draw[i,"relative_start"]),
                        as.numeric(features2draw[i,"relative_end"]),
                        as.numeric(features2draw[i,"relative_end"]))
          polygon_y = c(ymid,ymin,ymax)
        }
        #print(data.frame(polygon_x,polygon_y))
        return(data.frame(polygon_x,polygon_y))
      }))
    }
    #print(o)
    return(o)
  })

  if (!is.na(thisstart)) {
    thisdf$count = thisdf$count %>% filter(relative_pos > thisstart)
    # trim rect df
    # detect overlaps
    rectdf = rectdf %>% filter(xmax > thisstart)
    rectdf = rectdf %>% mutate(xmin = ifelse(xmin<thisstart,thisstart,xmin))
  }
  if (!is.na(thisend)) {
    thisdf$count = thisdf$count %>% filter(relative_pos < thisend)
    # trim rect df
    # detect overlaps
    rectdf = rectdf %>% filter(xmin < thisend)
    rectdf = rectdf %>% mutate(xmax = ifelse(xmax>thisend,thisend,xmax))
  }
  p <- ggplot()
  if (useReps) {
    for (i in unique(thisdf$count$rep)) {
        thisdata = thisdf$count %>% filter(rep==i)
        p <- p + geom_bar(data=thisdata,
                         aes(x = relative_pos, y = count, fill=segment),
                         stat="identity", position = "stack",
                         alpha=0.75) +
                 facet_grid(strain~targetSegment)
                            #,ncol=2,dir = "h")
    }
  } else {

  }

  p <- p + xlab("Position") + ylab("Count") +
        theme_few()  +
        geom_rect(data=rectdf,
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                    fill=segment)) +
        scale_fill_gradientn(colours = rev(rainbow(dim(segmentTable)[1],start=6/6,end=5/6)),
                        limits=c(0,dim(segmentTable)[1]),name = "Segment")
  if (logy) {
    p = p + scale_y_log10()
  }
  if (plot_genes) {
    # add gene annotations
    for (i in seq(1,dim(polygon_dfs)[1], by=3)) {
      p <- p + geom_polygon(data=polygon_dfs[i:(i+2),],
                          aes(x = polygon_x, y=polygon_y),fill="white")
    }
  }
  p
}
