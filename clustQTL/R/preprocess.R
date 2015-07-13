#' Converts GRanges object to matrix with real coordinates
#'
#' Some more info
#'
#' @param x An \code{m x n} matrix of peak intensities where the rows, \code{m},
#'  contains individuals, (e.g., strains)
#'  and the columns, \code{n}, contains observations (e.g., positions)
#' @param nresamples Number of resamples to compute.
#' @param alpha Significance level at which to reject the null hypothesis, Default = 0.01
#' @return A count value at which the null hypothesis can be rejected given \code{alpha}.
#' @examples
#' # S. cerevisiae chr01 3' isoform counts
#' x = subsetByOverlaps(tx3_counts_chr01,tx_3utr_annotations_yeast_chr01[1])
#' x = t(as.matrix(mcols(x)))
#' find_sigCounts(x)
#' @export
#'
granges2matrix = function(x, annot = NULL) {
  starts = start(x)
  if (!is.null(annot)) {
    m = matrix(0,nrow = length(seq(start(annot),end(annot))) ,
                               ncol = length(colnames(mcols(x))),
                               dimnames = list(seq(start(annot),end(annot)),
                                               colnames(mcols(x))))
  } else {
    m = matrix(0,nrow = length(seq(min(starts),max(starts))) ,ncol = length(colnames(mcols(x))),
               dimnames = list(seq(min(starts),max(starts)),colnames(mcols(x))))
  }
  tmp_m = as.data.frame(mcols(x))
  tmp_m = as.matrix(tmp_m)
  rownames(tmp_m) = starts
  m[rownames(tmp_m),colnames(tmp_m)] = tmp_m[rownames(tmp_m),colnames(tmp_m)]
  return(m)
}

#' Find significant counts in data matrix by resampling
#'
#' The function \code{\link{find_sigCounts}} returns a count cutoff at a specified
#'  significance value, \code{alpha}, by resampling. The null model specifies a uniform random
#'  distribution of counts across the sampling interval. The returned value signifies
#'  the count value at which the null hypothesis can be rejected given some confidence
#'  value, \code{alpha}
#'
#' @param x An \code{m x n} matrix of peak intensities where the rows, \code{m},
#'  contains individuals, (e.g., strains)
#'  and the columns, \code{n}, contains observations (e.g., positions)
#' @param nresamples Number of resamples to compute.
#' @param alpha Significance level at which to reject the null hypothesis, Default = 0.01
#' @return A count value at which the null hypothesis can be rejected given \code{alpha}.
#' @examples
#' # S. cerevisiae chr01 3' isoform counts
#' x = subsetByOverlaps(tx3_counts_chr01,tx_3utr_annotations_yeast_chr01[1])
#' x = t(as.matrix(mcols(x)))
#' find_sigCounts(x)
#' @export
#'
find_sigCounts = function(x, nresamples = 1000, alpha = 0.01) {
  countSamples = unlist(lapply(1:nresamples,function(i){
    colSums(do.call(rbind,lapply(1:dim(x)[1],function(i){sample(x[i,])})))
  }))
  quantile(countSamples,1-alpha)
}
