#' Permute
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return Modified mrk
#' @examples
#' format4manhattan()
#' @export
#'
format4manhattan = function( qtls, mrk, nresample = 10000 ) {
  mrk$p = -log10( qtls[,"pval"]+ 1/nresample )
  return(mrk)
}

#' Permute
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return plot
#' @examples
#' plotManhattan()
#' @export
#'
plotManhattan = function( qtls, mrk, qtl_name = "", trx_annot = NULL,
                          cutoff = 3, gene_annot_range = c(5000,5000),... ) {
  library(ggbio)
  mrk2 = format4manhattan( qtls,mrk )
  if ( (qtl_name!="") & ( !is.null(trx_annot) ) ) {
    # annotate gene
    trx_info = trx_annot[ which(trx_annot$Name == qtl_name), ]
    trx_granges = GRanges(seqnames=seqnames(trx_info),
                          ranges=IRanges(start(ranges(trx_info))-gene_annot_range[1],
                                         end(ranges(trx_info))+gene_annot_range[2]))
    names(trx_granges) = qtl_name
    plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                    ylab="-log10(pval)",main=qtl_name,ylim=c(0,4.5),
                    highlight.gr = trx_granges,...)
  } else {
    plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                    ylab="-log10(pval)",main=qtl_name,ylim=c(0,4.5),...)
  }

  #return(mrk2)
}
