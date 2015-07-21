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
                          cutoff = 3, gene_annot_range = c(5000,5000), chr = NULL,... ) {
  library(ggbio)
  if (is.null(rownames(qtls))) {
    # assume they are in correct order
    rownames(qtls) = names(mrk)
  }
  mrk2 = biovizBase::transformToGenome(format4manhattan( qtls,mrk ),0)
  if ( (qtl_name!="") & ( !is.null(trx_annot) ) ) {
    # annotate gene
    trx_info = trx_annot[ which(trx_annot$Name == qtl_name), ]
    trx_granges = GRanges(seqnames=seqnames(trx_info),
                          ranges=IRanges(start(ranges(trx_info))-gene_annot_range[1],
                                         end(ranges(trx_info))+gene_annot_range[2]))
    names(trx_granges) = qtl_name
    p = plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                    ylab="-log10(pval)",main=qtl_name,ylim=c(0,4.5),
                    highlight.gr = trx_granges,...)
  } else {
   p = plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                    ylab="-log10(pval)",main=qtl_name,ylim=c(0,4.5),...)

  }
  if (!is.null(chr)) {
    # chr functionality still not working correctly
    # some problems with chr breaks
    # only plot a specific chromosme
    l1 = min(mcols(mrk2[seqnames(mrk)==chr])$.start)
    l2 = max(mcols(mrk2[seqnames(mrk)==chr])$.end)
    #print(paste(l1,l2))
    p = p + xlim(l1,l2)
  }
  p = p + theme(axis.text.x=element_text(angle=-45, hjust=0))
  p
  return(p)
}
