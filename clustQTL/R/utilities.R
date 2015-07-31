#' Find QTL peaks
#'
#' The function \code{\link{findQTLPeaks}}
#'
#' @param qtls
#' @return Data Frame of peaks
#' @export
#'
findQTLPeaks = function(qtls, mrk, pcutoff = .05, peak_sigma = 25, peak_threshold=1,...) {
  # find peaks
  library(Peaks)
  library.dynam('Peaks', 'Peaks', lib.loc=NULL)
  p = as.numeric(qtls[,"pval"])
  locs = c(which(p>pcutoff),which(is.na(p)))
  p = -log10(p)
  p[locs] = 0
  names(p) = names(mrk)
  mrk2 = mrk
  mrk2$p = p
  peaks = SpectrumSearch(mcols(mrk2)$p,sigma=peak_sigma,threshold=peak_threshold)
  # peaks does not actually find QTL peak
  # to find this look for the QTL in a region
  # of peak with highest -log10 pval
  # this is actual "QTL" location
  actual_peaks = unlist(lapply(peaks$pos,function(i){
    i1 = i-peak_sigma*2
    if (i1<1) {
      i1 = 1
    }
    i2 = i+peak_sigma*2
    if (i2>length(p)) {
      i2 = length(p)
    }
    qtls_range = p[i1:i2]
    o = names(which(qtls_range==max(qtls_range)))
    o = o[median(seq(1,length(o)))]
    return(o)
  }))
  to_r = mrk2[actual_peaks]
  return(to_r)
}

#' getGenotype
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return genotypes
#' @export
#'
getGenotypes = function(marker,geno) {
  geno[marker,]
}

#' Cosine similarity
#'
#' Calculate cosine similarity between two vectors, \code{x} and \code{y}
#'
#' @param x A vector
#' @param y Another vector
#' @return Cosine similarity
#' @examples
#' x = rnorm(100,1)
#' y = rnorm(100,1)
#' cosine(x,y)
#' @export
#'
cosine = function(x,y) {
  x %*% y / sqrt(x%*%x * y%*%y)
}

#' cosineDist
#'
#' Calculate the cosine similarity between all rows of data matrix
#'  \code{x}
#'
#' @param x A matrix
#' @param na.rm Remove rows that generate NaN values
#' @return 1-cosine similarity matrix
#' @examples
#' x = do.call(rbind,lapply(1:5,function(i){rnorm(5,1)}))
#' cosineDist(x)
#' @export
cosineDist = function(x, na.rm=TRUE) {
  to_r = do.call(rbind,lapply(seq(1,dim(x)[1]),function(i){
    sapply(seq(1,dim(x)[1]),function(j){
      1-cosine(x[i,],x[j,])
    })
  }))
  # maintain names
  rownames(to_r) = rownames(x)
  colnames(to_r) = rownames(x)
  if (na.rm) {
    # clean up the cosineDist matrix by removing
    # find rows with all na values
    to_remove = which(apply(to_r,1,function(i)sum(is.nan(i))==length(i)))
    if (length(to_remove)>0) {
      to_r = to_r[-to_remove,-to_remove]
    }
    return(to_r)
  }
}

#' combineMarkers
#'
#' Combine markers if available genotypes at adjacent locations
#' are identical
#'
#' @param x A matrix
#' @param na.rm Remove rows that generate NaN values
#' @return 1-cosine similarity matrix
#' @export
combineMarkers = function(genotypes, markers, limit_strains = NULL, limit_markers = NULL,
                          na.rm = TRUE, rm_type = c("marker","strain")[1], collpase_markers = TRUE,
                          marker_rename = FALSE, impute_markers = TRUE, clean_markers = TRUE) {
  assertthat::assert_that(class(markers) == "GRanges" | class(markers) == "data.frame")
  if (class(markers) == "data.frame") {
    markers = GenomicRanges::makeGRangesFromDataFrame(markers,
       keep.extra.columns=TRUE,
       ignore.strand=TRUE)
  }
  if (is.null(names(markers))) {
    names(markers) = paste("mrk", seq(1,length(markers)), sep="_")
  }
  if (is.null(rownames(genotypes))) {
    rownames(genotypes) = names(markers)
  }
  genotypes = genotypes[order(as.numeric(gsub("mrk_", "", names(markers)))),]
  assertthat::assert_that(dim(genotypes)[1] == length(markers))

  if (clean_markers) {
    # clean markers according to Bloom et al 2013
    # marker must be called in 99% of segregants
    cat("...Removing markers genotyped in <95% of strains\n")
    tor = apply(genotypes,1,function(i){1-sum(is.na(i))/length(i)})
    tor = which(tor < 0.95)
    if (length(tor)>0) {
      genotypes = genotypes[-tor,]
    }
    # allele frequence should not be greater than 55% or less than 45%
    #tor = apply(genotypes,1,function(i){1-sum(i==1, na.rm=T)/length(i[!is.na(i)])})
    #tor = c(which(tor < 0.45), which(tor > 0.55))
    # replcaed 30/07/2015 with binomial test
    cat("...Removing markers with imbalanced allele frequency\n")
    tor = apply(genotypes,1,function(i){
      if (sum(!is.na(i)) == 0) {
        return(0)
      } else {
        binom.test(sum(i==1,na.rm=T),sum(!is.na(i)))$p.value
      }
    })
    tor = p.adjust(tor,method="BH")
    tor = which(tor<=0.05)
    if (length(tor)>0) {
      genotypes = genotypes[-tor,]
    }
    markers = markers[rownames(genotypes)]
    # for yeast genome genotypes_S288c_R64.rda
    # this leaves some markers with as many as 10 NA values
  }

  # record suspicious strains for user
  tor = apply(genotypes,2,function(i){1-sum(is.na(i))/length(i)})
  tor = which(tor<0.95)
  bad_strains = names(tor)

  if (impute_markers) {
    cat("...Imputing missing genotypes\n")
    #Fake phenotype data to make cross object
    phenotype = cbind(rep(0,dim(genotypes)[2]), rep(0,dim(genotypes)[2]))
    dimnames(phenotype) = list(colnames(genotypes),c("1","2"))
    genotype_subset = t( genotypes[ ,rownames( phenotype ) ] )
    # add chr num to markers
    genotype_subset = rbind( gsub('chr','', as.character( GenomicRanges::seqnames( markers[ colnames( genotype_subset ) ] ) ) ),
                             genotype_subset )
    genotype_subset = cbind( c( NULL, rownames( genotype_subset ) ), genotype_subset )
    colnames( genotype_subset )[ 1 ] = 'id'
    # add id column to phen data
    phenotype = cbind( phenotype, rownames( phenotype )  )
    colnames( phenotype )[ dim( phenotype )[2] ] = 'id'
    write.table( genotype_subset, file =  ".tmpgen", sep = ",", col.names = T, row.names = F, quote=F )
    write.table( phenotype, file =  ".tmpphen", sep = ",", col.names = T, row.names = F, quote=F )
    # read.cross to import data into rqtl
    cross = try(qtl::read.cross(format = "csvs", ".", genfile = ".tmpgen" ,
                                  phefile=  ".tmpphen", genotypes = c("1","2"),estimate.map=F))
    if ( class(cross)[1]!="try-error" ) {
      # clean up
      file.remove(c(".tmpgen",".tmpphen"))
      # impute by verterbi
      cross = qtl::fill.geno(cross,method="argmax")
      new_genotype = qtl::pull.geno(cross)
      rownames(new_genotype) = colnames(genotypes)
      genotypes = t(new_genotype)
    } else {
      cat("WARNING: Could not impute genotype \n")
    }
  }

  if (!is.null(limit_strains)) {
    assertthat::assert_that(sum(limit_strains%in%colnames(genotypes))>0)
    cat("...Limiting genotypes by user-supplied strains\n")
    if (sum(!limit_strains%in%colnames(genotypes))>0) {
      cat("WARNING: Not all limit_strains are in colnames(genotypes)")
    }
    genotypes = genotypes[, intersect(limit_strains, colnames(genotypes))]
  }
  if (!is.null(limit_markers)) {
    assertthat::assert_that(sum(limit_markers%in%names(markers))>0)
    cat("...Limiting genotypes by user-supplied markers\n")
    if (sum(!limit_markers%in%colnames(genotypes))>0) {
      cat("WARNING: Not all limit_markers are in names(markers)")
    }
    markers = markers[intersect(limit_markers, colnames(markers))]
    genotypes = genotypes[intersect(limit_markers, colnames(markers)),]
  }

  if(na.rm == TRUE) {
    # remove any strain and marker with NA values
    if (rm_type == "marker") {
      tor = apply(genotypes,1,function(i)sum(is.na(i)))
      if (sum(tor>0)>0) {
        genotypes = genotypes[tor>0,]
        markers = markers[tor>0]
      }
    } else if (rm_type == "strain") {
      tor = apply(genotypes,2,function(i)sum(is.na(i)))
      if (sum(tor>0)>0) {
        genotypes = genotypes[tor>0,]
        markers = markers[tor>0]
      }
    }
  }

  if(collpase_markers) {
    # collapse markers if genotypes are all the same
    # collapsed markers will take name of first marker by default
    # dim(genotypes)[1]
    cat("...Collapsing markers with identical genotypes. This may take some time...\n")
    new_markers = GenomicRanges::GRanges()
    i = dim(genotypes)[1]-1000
    pb <- txtProgressBar(min = 1, max =   dim(genotypes)[1], style = 3)
    while(i <=  dim(genotypes)[1]) {
      setTxtProgressBar(pb, i)
      i1 = i + 1
      if (i1 <= dim(genotypes)[1]) {
        while(
          # same genotype
          sum(diff(rbind(genotypes[i,],genotypes[i1,]))==0)==length(genotypes[i,]) &
          # same chromosome
          as.logical(GenomicRanges::seqnames(markers[i])==GenomicRanges::seqnames(markers[i1])) &
          # not out of bounds
          i1 <= dim(genotypes)[1]
        ) {
          i1 = i1 + 1
          #print(i1)
        }
        new_marker = markers[i,]
        GenomicRanges::end(GenomicRanges::ranges(new_marker)) =
          GenomicRanges::end(GenomicRanges::ranges(markers[i1-1,]))
        new_markers = c(new_markers,new_marker)
      } else {
        new_marker = markers[i,]
        new_markers = c(new_markers,new_marker)
      }
      # print(paste(i,i1))
      # adjust marker
      i = i1
    }
    close(pb)
    GenomicRanges::mcols(new_markers) = NULL
    markers = new_markers
    genotypes = genotypes[names(markers),]
  }

  if(marker_rename == TRUE) {
    names(markers) = paste("mrk", seq(1,length(markers)), sep="_")
    rownames(genotypes) = names(markers)
  }
  return(list(genotypes = genotypes, markers = markers, suspect_strains = bad_strains))
}

#tmp = combineMarkers(geno,mrk,limit_strains=colnames(genotypes)[1:10])
