#' Format clustQTL results for plotting
#'
#' This function will format results from running clustQTL for plotting
#'  as a Manhattan plot with \code{\link{plotManhattan}}. Plotting
#'  is done with plotGradLinear from the ggbio package. This function is
#'  not usually called directly. It is called by \code{\link{plotManhattan}}
#'
#' @param qtls A two column matrix with each row containing a pvalue for
#'  every marker in \code{mrk}
#' @param mrk A Granges object containing the location of every genetic marker
#'  tested by clustQTL
#' @return GRanges object containing an extra metacolumn "p"
#'  containing the -log10(pvalue) at each marker. This object is suitable
#'  for plotting with \code{\link{plotManhattan}}
#' @examples
#' format4manhattan()
#' @export
#'
format4manhattan = function(qtls, mrk) {
  mrk$p = -log10(qtls[,"pval"])
  return(mrk)
}

#' Plot clustQTL results
#'
#' The function \code{\link{cluster}}
#'
#' @param qtls A two column matrix with each row containing a pvalue for
#'  every marker in \code{mrk}
#' @param mrk A Granges object containing the location of every genetic marker
#'  tested by clustQTL
#' @param main String used for plot title
#' @param trx_annot GRanges. \code{main} show
#' @param cutoff -log10(pval) signficance cutoff value for QTLs. Used to draw
#'  horizontal line across plot
#' @param gene_annot_range
#' @param cutoff
#' @return plot
#' @examples
#' plotManhattan()
#' @export
#'
plotManhattan = function( qtls, mrk, main = "", trx_annot = NULL,
                          cutoff = 3, gene_annot_range = c(5000,5000), chr = NULL,... ) {
  library(ggbio)
  if (is.null(rownames(qtls))) {
    # assume they are in correct order
    rownames(qtls) = names(mrk)
  }
  mrk2 = biovizBase::transformToGenome(format4manhattan( qtls,mrk ),0)
  if ( (qtl_name!="") & ( !is.null(trx_annot) ) ) {
    # annotate gene
    trx_info = trx_annot[ which(trx_annot$Name == main), ]
    trx_granges = GRanges(seqnames=seqnames(trx_info),
                          ranges=IRanges(start(ranges(trx_info))-gene_annot_range[1],
                                         end(ranges(trx_info))+gene_annot_range[2]))
    names(trx_granges) = main
    p = plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                    ylab="-log10(pval)",main=main,ylim=c(0,4.5),
                    highlight.gr = trx_granges,...)
  } else {
   p = plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                    ylab="-log10(pval)",main=maine,ylim=c(0,4.5),...)

  }
  if (!is.null(chr)) {
    # chr functionality still not working correctly
    # some problems with chr breaks so it plots incorrect region
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
