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
findQTLPeaks = function( qtls, mrk, cutoff = 3 ) {
  # find top two peaks
  mrk2 = format4manhattan( qtls, mrk )
  candidates = which( mrk2$p>cutoff)
  ind = 1
  final_candidates = cbind(candidates[ind],mrk2$p[candidates[ind]])
  ind = ind + 1
  while (ind<=length(candidates)) {
    if (diff(candidates)[ind-1]<10) {
      # single peak
      # only replace val if higher
      if ( mrk2$p[candidates[ind]] > final_candidates[ dim(final_candidates)[1],2] ) {
        final_candidates[ dim(final_candidates)[1],] = cbind(candidates[ind],mrk2$p[candidates[ind]])
      }
    } else {
      final_candidates = rbind(final_candidates,cbind(candidates[ind],mrk2$p[candidates[ind]]))
    }
    ind = ind + 1
  }
  to_r = mrk2[paste("mrk",final_candidates[,1],sep="_"),]
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

