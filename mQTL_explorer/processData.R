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
.local = FALSE
if (system("hostname",intern=T) == "mac-steinmetz55.embl.de") {
  print("yes")
  .local = TRUE
} else {
  print(system("hostname"))
}

# Import packages ---------------------------------------------------
library(funqtl)
library(GenomicRanges)
library(rtracklayer)
library(jsonlite)

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
if (.local) {
  DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"
  WDIR = "/Users/brooks/Sites/JBrowse-1.11.6_mQTL"
  EDIR = "/Users/brooks/Documents/steinmetz_local/genphen/transcriptome"
} else {
  DDIR = "/g/steinmetz/brooks/genphen/metabolome/qtls"
  WDIR = "/var/www2/html/mQTL"
  EDIR = "/g/steinmetz/brooks/genphen/transcriptome"
}


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
  trackList = fromJSON(file.path(WDIR, "data/trackList.json"))

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
  toadd = lapply(names(data), function(i) {
    if (!is.null(data[[i]]$qtl)) {
      g$score = unlist(lapply(g$revmap,function(j){mean(data[[i]]$qtl[j,4], na.rm=T)}))
      export.bw(object = g,
                con = file.path(WDIR,paste("data/tracks/mqtls/",i,".bw",sep="")))
      # add trackList entry
      if (!(i %in% trackList$tracks[,"key"])) {
        newRow = fromJSON(sprintf(MQTL, i, i, i))
        return(newRow)
      } else {
        return(newRow)
      }
    return(newRow)
    }
  })
  toadd = toadd[!sapply(toadd,is.null)]
  if (length(toadd)>1) {
    trackList$tracks = rbind.pages(list(trackList$tracks,fromJSON(toJSON(toadd,pretty=T,auto_unbox=T))))
  } else if (length(toadd)>1) {
    cat("ERROR:Cannot add only one track.\n")
  }
  write(toJSON(trackList, pretty= TRUE, auto_unbox = T), file = file.path(WDIR, "data/trackList.json"))
}

  if (.makeeQTL) {
    # make eQTL tracks
    cat("Writing eQTL BigWig files and adding to JBrowse trackList config file. This may take some time.\n")
    load(file.path(EDIR,"eQTL.rda"))
    seqinfo = Seqinfo(genome="sacCer3")
    trackList = fromJSON(file.path(WDIR, "data/trackList.json"))

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
    pb = txtProgressBar(min = 3, max = dim(eQTL$qtls)[2], style = 3)
    toadd = lapply(3:dim(eQTL$qtls)[2], function(i) {
      #print(i)
      #setTxtProgressBar(pb, i)
      if (sum(is.na(eQTL$qtls[,i])) == 0) {
        i = as.numeric(i)
        ii = colnames(eQTL$qtls)[i]
        ii = sub("\\.2C","-",ii)
        ii = sub("\\.","-",ii)
        tof = file.path(WDIR,paste("data/tracks/eqtls/",ii,".bw",sep=""))
        if (!file.exists(tof)) {
          print(ii)
          g$score = unlist(lapply(g$revmap,function(j){mean(eQTL$qtls[j,i], na.rm=T)}))
          export.bw(object = g,
                    con = tof)
        }
        # add trackList entry
        if (!(ii %in% trackList$tracks[,"key"])) {
          newRow = fromJSON(sprintf(EQTL, ii, ii, ii))
          return(newRow)
        } else {
          return(NULL)
        }
      }
    })
    toadd = toadd[!sapply(toadd,is.null)]
    close(pb)
    if (length(toadd)>1) {
      trackList$tracks = rbind.pages(list(trackList$tracks,fromJSON(toJSON(toadd,pretty=T,auto_unbox=T))))
    } else if (length(toadd)>1) {
      cat("ERROR:Cannot add only one track.\n")
    }
    write(toJSON(trackList, pretty= TRUE, auto_unbox = T), file = file.path(WDIR, "data/trackList.json"))
}
