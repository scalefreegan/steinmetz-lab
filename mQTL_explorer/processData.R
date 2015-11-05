#-------------------------------------------------------------------#
# Shiny interface for ploting/exploring mQTLs
# 
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "WTFPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.makeBW = FALSE
.makeeQTL = FALSE
# Import packages ---------------------------------------------------
library(funqtl)
library(GenomicRanges)
library(rtracklayer)
library(RJSONIO)

# Templates ---------------------------------------------------
# mqtl
MQTL = strwrap('{"category": "mQTL",
  "label"         : "%s",
  "key"           : "%s",
  "storeClass"    : "JBrowse/Store/SeqFeature/BigWig",
  "urlTemplate"   : "tracks/mqtls/%s.bw",
  "type"          : "JBrowse/View/Track/Wiggle/Density",
  "bicolor_pivot" : "zero",
  "min_score"     : "0",
  "max_score"     : "10",
  "style": {
    "pos_color": "purple",
    "neg_color": "green"
  }
}')

EQTL = strwrap('{"category": "eQTL",
               "label"         : "%s",
               "key"           : "%s",
               "storeClass"    : "JBrowse/Store/SeqFeature/BigWig",
               "urlTemplate"   : "tracks/eqtls/%s.bw",
               "type"          : "JBrowse/View/Track/Wiggle/Density",
               "bicolor_pivot" : "zero",
               "min_score"     : "0",
               "max_score"     : "10",
               "style": {
               "pos_color": "green",
               "neg_color": "purple"
               }
               }')

# tmp resource location / will be changed
DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"
WDIR = "/Users/brooks/Sites/JBrowse-1.11.6_mQTL"
EDIR = "/Users/brooks/Documents/steinmetz_local/genphen/transcriptome"

load(file.path(DDIR,"endometabolite_full_12102015.rda"))
load(file.path(DDIR,"genotypes_S288c_R64.rda"))
load(file.path(DDIR,"mQTLs_comball_funqtl_2014.rda"))

data = lapply(seq(1,length(mQTLs_funqtl_2014)), function(i){
  if (class(mQTLs_funqtl_2014[[i]])=="try-error") {
    return(NULL)
  } else {
    #print(i)
    o = list()
    o[["qtl"]] = mQTLs_funqtl_2014[[i]][["qtls_alt"]]
    o[["permout"]] = mQTLs_funqtl_2014[[i]]$permout
    return(o)
  }
})
names(data) = names(mQTLs_funqtl_2014)

# remove NULLs
data = data[!sapply(data,is.null)]

if (.makeBW) {
  seqinfo = Seqinfo(genome="sacCer3")
  trackList = fromJSON(file.path(WDIR, "data/trackList.json"),flatten = TRUE)
  
  # make mQTL tracks
  nn = sapply(as.character(levels(seqnames(mrk))),function(i){
    paste(substr(i,1,3),as.roman(substr(i,4,5)),sep="")
  })
  mrk_mod = renameSeqlevels(mrk,nn)
  mrk_mod = keepSeqlevels(mrk_mod,unique(nn))
  mrk_mod$score = 1
  g = GenomicRanges::GRanges(
    seqnames = seqnames(mrk_mod),
    ranges = ranges(mrk_mod),
    strand = strand(mrk_mod),
    seqinfo = seqinfo
  )
  g = reduce(g,with.revmap=T)
  toadd = data.frame()
  for (i in names(data)) {
    if (!is.null(data[[i]]$qtl)) {
      g$score = unlist(lapply(g$revmap,function(j){mean(data[[i]]$qtl[j,4], na.rm=T)}))
      export.bw(object = g, 
                con = file.path(WDIR,paste("data/tracks/mqtls/",i,".bw",sep="")))
      # add trackList entry
      if (!(i %in% trackList$tracks[,"key"])) {
        newRow = as.data.frame(fromJSON(sprintf(MQTL, i, i, i)))
        rownames(newRow) = i
        toadd = rbind(toadd,newRow)
      }
    }
  }
  if (length(toadd)>0) {
    trackList$tracks = merge(trackList$tracks, toadd, all = TRUE, sort = FALSE)
  }
  
  if (.makeeQTL) {
    # make eQTL tracks
    cat("Writing eQTL BigWig files and adding to JBrowse trackList config file. This may take some time.\n")
    load(file.path(EDIR,"eQTL.rda"))
    toadd = data.frame()
    pb = txtProgressBar(min = 0, max = dim(eQTL$qtls)[2]-2, style = 3)
    #dim(eQTL$qtls)[2]
    for (i in 3:10) {
      setTxtProgressBar(pb, i)
      if (sum(is.na(eQTL$qtls[,i])) == 0) {
        g$score = unlist(lapply(g$revmap,function(j){mean(eQTL$qtls[j,i], na.rm=T)}))
        i = colnames(eQTL$qtls)[i]
        export.bw(object = g, 
                  con = file.path(WDIR,paste("data/tracks/eqtls/",i,".bw",sep="")))
        # add trackList entry
        if (!(i %in% trackList$tracks[,"key"])) {
          # get a template row
          newRow = trackList$tracks[1,] 
          newRow[] = NA
          newRow = as.data.frame(fromJSON(sprintf(EQTL, i, i, i)))
          rownames(newRow) = i
          toadd = rbind(toadd,newRow)
        }
      }
    }
    for (x in setdiff(colnames(trackList$tracks),colnames(toadd))) {
      print(x)
      if (!is.null(dim(trackList$tracks[[x]])[2]) ) {
        toadd[[x]] = matrix(NA, nrow = dim(toadd)[1], ncol = dim(trackList$tracks[[x]])[2])
        colnames(toadd[[x]]) = colnames(trackList$tracks[[x]])
      } else {
        toadd[[x]] = NA
      }
    }
    trackList$tracks = rbind(trackList$tracks, toadd)
  }
  
  # write json trackList
  trackList = lapply(trackList,as.data.frame)
  write(toJSON(trackList, pretty= TRUE), file = file.path(WDIR, "data/trackList.json"))
}
