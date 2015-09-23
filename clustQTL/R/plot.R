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
#' @export
#'
format4manhattan = function(qtls, mrk, qqman = TRUE, midpoint = TRUE, restrict = NULL, ...) {
  if (qqman) {
    if (midpoint) {
      # take midpoint instead of entire range
      bp =  apply(cbind(GenomicRanges::start(GenomicRanges::ranges(mrk)),
                       GenomicRanges::end(GenomicRanges::ranges(mrk))),1,function(i){round(mean(i))})
      snp = names(mrk)
      chr = as.numeric(gsub("chr","",GenomicRanges::seqnames(mrk)))
      p = qtls[,"pval"]
    } else {
      # Not implemented
      #bp =  cbind(GenomicRanges::start(GenomicRanges::ranges(mrk)),
      #          GenomicRanges::end(GenomicRanges::ranges(mrk)))
      #snp = rep(names(mrk),length(bp))
      #chr = rep(GenomicRanges::seqnames(mrk),length(bp))
      #p = rep(qtls[,"pval"],length(bp))
    }
    to_r = data.frame(SNP = snp,CHR = chr, BP = bp, P = p)
    if (!is.null(restrict)) {
      to_r = to_r %>% filter(.,SNP == restrict)
    }
    return(to_r)
  } else {
    mrk$p = -log10(qtls[,"pval"])
    return(mrk)
  }
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
#' @export
#'
plotManhattan = function( qtls, mrk, gene = "", trx_annot = NULL,
                          cutoff = 3, gene_annot_range = c(1000,1000), qqman = T, show = T, ... ) {
  if (qqman) {
    if (is.null(rownames(qtls))) {
      # assume they are in correct order
      rownames(qtls) = names(mrk)
    }
    mrk2 = format4manhattan( qtls, mrk, qqman = T)
    if (!is.null(trx_annot)) {
      trx_info = trx_annot[ which(trx_annot$Name == gene), ]
      trx_granges = GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(trx_info),
                                           ranges=IRanges::IRanges(
                                             GenomicRanges::start(GenomicRanges::ranges(trx_info))-gene_annot_range[1],
                                             GenomicRanges::end(GenomicRanges::ranges(trx_info))+gene_annot_range[2]))
      names(trx_granges) = gene
      gene_snps = mrk2 %>% filter(., CHR == as.numeric(gsub("chr","",GenomicRanges::seqnames(trx_granges))),
                                  BP >= GenomicRanges::start(GenomicRanges::ranges(trx_granges)),
                                  BP <= GenomicRanges::end(GenomicRanges::ranges(trx_granges)))
      to_r = qqman::manhattan(mrk2, highlight = gene_snps$SNP, show = show, ...)
    } else {
      to_r = qqman::manhattan(mrk2, show = show, ...)
    }
    return(to_r)
  } else {
    library(ggplot2)
    if (is.null(rownames(qtls))) {
      # assume they are in correct order
      rownames(qtls) = names(mrk)
    }
    mrk2 = format4manhattan( qtls, mrk, qqman = F )
    if ( (gene != "") & ( !is.null(trx_annot) ) ) {
      # annotate gene
      trx_info = trx_annot[ which(trx_annot$Name == gene), ]
      trx_granges = GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(trx_info),
                                           ranges=IRanges::IRanges(
                                             GenomicRanges::start(GenomicRanges::ranges(trx_info))-gene_annot_range[1],
                                             GenomicRanges::end(GenomicRanges::ranges(trx_info))+gene_annot_range[2]))
      names(trx_granges) = gene
      p = ggbio::plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                                 ylab="-log10(pval)",
                                 highlight.gr = trx_granges,...)
    } else {
      p = ggbio::plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,
                                 ylab="-log10(pval)",...)

    }
    p + theme(axis.text.x=element_text(angle=-45, hjust=0))
  }

}

#' Plot peak profiles
#'
#' The function \code{\link{cluster}}
#'
#' @param qtls A two column matrix with each row containing a pvalue for
#'  every marker in \code{mrk}
#' @return plot
#' @export
#'
plotPeakProfile = function(data, genotypes, marker, peak_sigma = 2, peak_threshold = 1) {
  library(reshape2)
  library(Peaks)
  library.dynam('Peaks', 'Peaks', lib.loc=NULL)
  library(gtable)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(ggplot2)
  data_mod = t(apply(data,1,function(i){SpectrumSearch(i,sigma=peak_sigma,threshold=peak_threshold)$y}))
  #data_mod = t(apply(data_mod,1,function(i){i/sum(i)}))
  colnames(data_mod) = colnames(data)
  data_long = melt(data, varnames=c("y","x"))
  data_long$genotype = genotypes[marker,][levels(data_long$y)[data_long$y]]
  data_long$panel = "B"
  colnames(data_long) = c("value","x","y","genotype","panel")
  g1 = data_long %>% filter(genotype==1)
  g2 = data_long %>% filter(genotype==2)
  s2n = c(seq(1,length(unique(levels(g1$value)[g1$value]))),seq(1,length(unique(levels(g1$value)[g2$value]))))
  names(s2n) = c(unique(levels(g1$value)[g1$value]),unique(levels(g2$value)[g2$value]))
  data_long$value = s2n[levels(data_long$value)[data_long$value]]
  data_long_mod = melt(data_mod, varnames=c("y","x"))
  data_long_mod$genotype = genotypes[marker,][levels(data_long_mod$y)[data_long_mod$y]]
  data_long_mod$panel = "A"
  data_long_mod$y = s2n[levels(data_long_mod$y)[data_long_mod$y]]
  d = rbind(data_long,data_long_mod)
  # remove 0 values
  # data_long = data_long[data_long$value!=0,]
  # data_long[data_long$value==0,]$value = NA

  #panelAylim = max(data_long_mod %>% filter(genotype==1) %>% group_by(genotype,x) %>% summarise(sum=sum(value)) %>% select(sum))
  p1 <- ggplot() +
    facet_grid(panel~.,scale="free_y",labeller=function(x,y){return("")}) +
    stat_summary(data=data_long_mod %>% filter(genotype==1), mapping=aes(y=value,x=x),geom="area",fun.y=sum) +
    #scale_y_continuous(limits=c(0, panelAylim)) +
   # geom_line(data=data_long_mod %>% filter(genotype==1), mapping=aes(y=value,x=x),stat="mean") +
    geom_tile(data = data_long %>% filter(genotype==1), mapping=aes(y = value,x = x,fill = y)) +
    scale_fill_gradient(low="black", high="white", guide = F, limits=c(0, max(data_long$y))) +
    labs(x = "Position", y = "Sum Counts") +
    ggtitle("Genotype 1")
  gt1 <- ggplot_gtable(ggplot_build(p1))
  gt1$heights[[3]] <- unit(.25, "null")
  gt1$grobs[[3]]$children[2]$axis$grobs[[1]]$label = ""
  gt1$grobs[[10]]$hjust = -1.5
  gt1$grobs[[10]]$vjust = -.25


  p2 <- ggplot() +
    facet_grid(panel~.,scale="free",labeller=function(x,y){return("")}) +
    stat_summary(data=data_long_mod %>% filter(genotype==2), mapping=aes(y=value,x=x),geom="area",fun.y=sum) +
    #scale_y_continuous(limits=c(0, panelAylim)) +
    #geom_line(data=data_long_mod %>% filter(genotype==2), mapping=aes(y=value,x=x)) +
    geom_tile(data = data_long %>% filter(genotype==2), mapping=aes(y=value,x=x,fill = y)) +
    scale_fill_gradient(low="black", high="white", limits=c(0, max(data_long$y)),name="Counts") +
    labs(x = "Position", y = "Sum Counts") +
    ggtitle("Genotype 2")
  gt2 <- ggplot_gtable(ggplot_build(p2))
  gt2$heights[[3]] <- unit(.25, "null")
  gt2$grobs[[3]]$children[2]$axis$grobs[[1]]$label = ""
  gt2$grobs[[10]]$hjust = -1.5
  gt2$grobs[[10]]$vjust = -.25

  grid.arrange(gt1,gt2,ncol=2)

}

#' Plot peak profiles
#'
#' The function \code{\link{cluster}}
#'
#' @param qtls A two column matrix with each row containing a pvalue for
#'  every marker in \code{mrk}
#' @return plot
#' @export
#'
writePeakProfile = function(data, genotypes, marker, peak_sigma = 2, peak_threshold = 1, 
                            f1 = "/Users/brooks/Sites/JBrowse-1.11.6/data/d1.bw", f2 = "/Users/brooks/Sites/JBrowse-1.11.6/data/d2.bw") {
  library(reshape2)
  library(Peaks)
  library.dynam('Peaks', 'Peaks', lib.loc=NULL)
  library(gtable)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(ggplot2)
  library(rtracklayer)
  mrk2geno = genotypes[marker,]
  #data_mod = as.matrix(mcols(data))
  #data_mod = t(apply(data,1,function(i){SpectrumSearch(i,sigma=peak_sigma,threshold=peak_threshold)$y}))
  strain_names = sapply(colnames(mcols(data)),function(i){
    i = strsplit(i,"_")[[1]][1]
    i = gsub("X","",i)
    if (nchar(i)==2) {
      i = paste("0",i,sep="")
    }
    if (nchar(i)==3) {
      return(i)
    } else {
      return(NULL)
    }
  })
  gen1 = which(strain_names%in%names(which(mrk2geno==1)))
  gen2 = which(strain_names%in%names(which(mrk2geno==2)))
  d1 = rowSums(as.matrix(mcols(data)[,gen1]))
  d2 = rowSums(as.matrix(mcols(data)[,gen1]))
  gr1 = data
  mcols(gr1) = NULL
  gr1$score = d1
  gr2 = data
  mcols(gr2) = d2
  export.bw(object = gr1, con = f1, compress=T)
}
