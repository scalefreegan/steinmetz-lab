#' Find QTL peaks
#'
#' The function \code{\link{findQTLPeaks}}
#'
#' @param qtls
#' @return Data Frame of peaks
#' @examples
#' findQTLPeaks()
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
#' @examples
#' getGenotypes()
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

