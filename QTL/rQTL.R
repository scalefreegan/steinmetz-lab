#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# Rqtl
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Alpha"

PCA <- function(mat) eigen(cov(apply(mat, 2, function(i) (i - mean(i))/sd(i))))

removePrincipalComponent <- function(
  # remove one or more principal components from data matrix
  matAdjust = 'Centered, variance scaled matrix',
  meanList = 'Column means of original (unadjusted) matrix',
  eigenList = 'List of eigenvalues and eigenvectors of adjust matrix covariance matrix',
  n = 'selected PC\'s to remove',
  specific_select = 'If True: n == 1:n, if False: just n\'th columns') {

  if (length(n) > ncol(matAdjust)) stop('N is higher than the number of PC\'s')
  if (!specific_select & length(n) > 1) stop('Use a single number when selecting up to n\'th PC')
  if (!specific_select) n <- 1:n

  to_r = t(eigenList$vectors[,-n] %*% (t(eigenList$vectors[,-n]) %*% t(matAdjust))) + (matrix(meanList, nrow = nrow(matAdjust), ncol = ncol(matAdjust), byrow=T))
  colnames(to_r) = colnames(matAdjust)
  return(to_r)
}

runQTL <- function(
    # Streamline rQTL
    # Format standard matrix/genotype data for rQTL
    # Run rQTL qtl::scanone analysis with optional methods,
    # permute_sig and pca
    genotype, # a n x m genotype matrix containing (eg 1,2) at 'n' genotype markers in 'm' strains
    phenotype, # a n x m phenotype matrix containing measurements of 'm' phenotypes in 'n' strains
    marker_info, # GRanges object with information about genetic markers (chromosome and location)
    permute = T, # compute significance of each QTL LOD by permutation
    pca = F, # maximize QTL detection by removing confounders using PCA
    permute_alpha = 0.05,
    permute.lod.cutoff = 2.5,
    save_file = "",
    return_cross = FALSE, # just return cross object,
    subset_genotype = F,
    method = "hk",
    n.cores = 24,
    n.perm = 1000,
    genotypes = c( "1","2" ),
    estimate.map = FALSE,
    ...
    ){
  # subset genotype data on metabolite data. transpose it.
  if ( !sum(colnames(genotype)%in%rownames(phenotype))>0 ) {
    cat("ERROR: Strain names in genotype matrix (columns) do not match strain names in phenotype matrix (rows\n")
    return(NULL)
  }
  if (sum(is.na(pheno))>0 && pca == TRUE) {
    cat("ERROR: PCA cannot handle missing/NA values\n")
    return(NULL)
  }
  if (subset_genotype) {
    genotype = t( genotype[ ,rownames( phenotype ) ] )
    # add chr num to markers
  } else {
    genotype = t(genotype)
  }
  genotype = rbind( gsub('chr','', as.character( GenomicRanges::seqnames( marker_info[ colnames( genotype ) ] ) ) ),
   genotype )
  genotype = cbind( c( NULL, rownames( genotype ) ), genotype )
  colnames( genotype )[ 1 ] = 'id'

  # add id column to phen data
  phenotype = cbind( phenotype, rownames( phenotype )  )
  colnames( phenotype )[ dim( phenotype )[2] ] = 'id'
  # write rqtl files
  write.table( genotype, file =  ".tmpgen", sep = ",", col.names = T, row.names = F, quote=F )
  write.table( phenotype, file =  ".tmpphen", sep = ",", col.names = T, row.names = F, quote=F )
  # read.cross to import data into rqtl
  cat("Creating rQTL cross object for genotype/phenotype data\n")
  genphen = try( qtl::read.cross( format = "csvs", ".", genfile = ".tmpgen" , phefile=  ".tmpphen", genotypes = genotypes, estimate.map = estimate.map, ... ) )
  if ( class(genphen)[1]!="try-error" ) {
    # clean up
    file.remove(c(".tmpgen",".tmpphen"))
  } else {
    cat("ERROR: Could not create rQTL cross object from genotype and phenotype data\nDo you have write permission in the current directory?\n")
    return(NULL)
  }
  genphen = qtl::calc.genoprob( genphen, step = 0 )
  if (return_cross) {
    # just return cross object, nothing else
    # to use in funQTL etc
    return(genphen)
  }
  if (pca) {
    cat("Removing principal components to increase number of detected QTLs\n")
    pc_removed = 0
    # as a first approximation, use lod score > 2.5 as a "true" QTL
    # since running permutation is expensive
    #phenotype = genphen$pheno[,colnames(genphen$pheno)!="id",drop=F]
    #tls = sum(unlist(parallel::mclapply( seq(1,dim(phenotype)[2]),function( i ) {sum(qtl::scanone(genphen, pheno.col = i, method = method)$lod >= 2.5)})))
    phes = which(colnames(genphen$pheno)!="id")
    qtls = qtl::scanone(genphen, pheno.col = phes, method = method)
    qtls = sum(qtls[,intersect(colnames(qtls),colnames(genphen$pheno))] >= permute.lod.cutoff)
    phenotype = genphen$pheno[,colnames(genphen$pheno)!="id",drop=F]
    # remove NAs
    phenotype = phenotype[!apply(phenotype,1,function(x){sum(is.na(x))>0}),]
    qtls_mod = 0
    pca <- PCA(phenotype)
    phenotype_mod <- removePrincipalComponent(
      matAdjust = apply(phenotype, 2, function(i) (i - mean(i))/sd(i)),
      meanList = apply(phenotype, 2, mean),
      eigenList = pca,
      n = 1,
      specific_select = TRUE
    )
    genphen_mod = genphen
    genphen_mod$pheno[rownames(phenotype_mod), colnames(genphen$pheno)!="id"] = phenotype_mod
    qtls_mod = qtl::scanone(genphen_mod, pheno.col = phes, method = method)
    qtls_mod = sum(qtls_mod[,intersect(colnames(qtls_mod),colnames(genphen_mod$pheno))] >= permute.lod.cutoff)
    while (qtls_mod>qtls) {
      pc_removed = pc_removed + 1
      phenotype = phenotype_mod
      qtls = qtls_mod
      pca <- PCA(phenotype)
      phenotype_mod <- removePrincipalComponent(
        matAdjust = apply(phenotype, 2, function(i) (i - mean(i))/sd(i)),
        meanList = apply(phenotype, 2, mean),
        eigenList = pca,
        n = 1,
        specific_select = TRUE
      )
      genphen_mod = genphen
      genphen_mod$pheno[rownames(phenotype_mod), colnames(genphen$pheno)!="id"] = phenotype_mod
      qtls_mod = qtl::scanone(genphen_mod, pheno.col = phes, method = method)
      qtls_mod = sum(qtls_mod[,intersect(colnames(qtls_mod),colnames(genphen_mod$pheno))] >= permute.lod.cutoff)
    }
    if (pc_removed>0) {
      pca <- PCA(phenotype)
      phenotype <- removePrincipalComponent(
        matAdjust = apply(phenotype, 2, function(i) (i - mean(i))/sd(i)),
        meanList = apply(phenotype, 2, mean),
        eigenList = pca,
        n = 1:pc_removed,
        specific_select = TRUE
      )
    }
    to_r = list()
    genphen$pheno[rownames(phenotype), colnames(genphen$pheno)!="id"] = phenotype
    phes = which(colnames(genphen$pheno)!="id")
    to_r$qtls = qtl::scanone(genphen, pheno.col = phes, method = method)
    to_r$pc_removed = pc_removed
    to_r$phenotype = phenotype
  } else {
    to_r = list()
    phes = which(colnames(genphen$pheno)!="id")
  	to_r$qtls = qtl::scanone(genphen, pheno.col = phes, method = method)
  }
  to_r$cross = genphen
  if (permute) {
    phes = which(colnames(genphen$pheno)!="id")
    # permuations to determine qtl sig cutoff
    to_r$qtls_permuted = qtl::scanone( cross, pheno.col = phes, n.perm = n.perm, n.cluster = n.cores)
  }
  if (save_file!="") {
    qtl_list = to_r
    save(qtl_list,file=save_file)
  }
  return(to_r)
}

plotlod_wlabs <- function(output, effects, y, y.text, ylab="Time", gap=25,
                    ncolors=251, horizontal=FALSE, ...)
{

    if(missing(y)) {
        y <- 1:(ncol(output)-2)
    }

    if(missing(y.text)) {
        y.text <- y
    }

    templod <- as.matrix(output[,-(1:2)])
    maxlod <- max(templod, na.rm=TRUE)
    zlim <- c(0, maxlod)
    val <- sqrt(seq(0, 1, len=ncolors))
    col <- rgb(1, rev(val), rev(val))
    if(!missing(effects)) {
        if(!all(dim(effects) == dim(templod)))
            stop("dim(effects) doesn't conform to dim(output)")
        templod[effects < 0] <- templod[effects<0] * -1
        zlim <- c(-maxlod, maxlod)
        ncolors=(ncolors-1)/2+1
        val <- sqrt(seq(0, 1, len=ncolors))
        col <- c(rgb(val, val, 1), rgb(1, rev(val), rev(val))[-1])
    }

    uchr <- unique(output[,1])
    pos <- NULL
    lod <- NULL
    chr <- vector("list", length(uchr))
    names(chr) <- uchr
    off.end <- 0.5 # spacing off the ends of the chromosomes
    for(i in uchr) {
        temppos <- output[output[,1]==i,2]
        temppos <- temppos - min(temppos)
        temppos <- rowMeans(cbind(c(temppos[1]-off.end, temppos),
                                  c(temppos, max(temppos)+off.end)))
        if(is.null(pos)) {
            pos <- temppos
            lod <- templod[output[,1]==i,]
        }
        else {
            temppos <- max(pos)+gap + temppos
            pos <- c(pos, temppos)
            lod <- rbind(lod, rep(0, ncol(lod)), templod[output[,1]==i,])
        }
        chr[[i]] <- c(min(temppos), max(temppos))
        chr[[i]] <- c(chr[[i]], mean(chr[[i]]))
    }




    if(horizontal == FALSE) {

        ## add more space between regend and plot
        pos <- c(pos, max(pos) + gap)
        lod <- rbind(lod, rep(0, ncol(lod)))

        ## plot
        image.plot(pos, y, lod, xaxt="n", ylab=ylab, xlab="",
              zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n", axes=F)
        image(pos, y, lod, axes=FALSE, add=TRUE, pos, y, lod, xaxt="n", ylab=ylab, xlab="",
              zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n", axes=F)
        mtext(text=y.text, side=2, line=0.3, at=y, las=1, cex=0.8)
        title(xlab="Chromosome", mgp=c(2, 0, 0))
        u <- par("usr")
        yd <- 0.04
        for(i in seq(along=chr)) {
            rect(chr[[i]][1]-.5, u[3], chr[[i]][2]+.5, u[3]-diff(u[3:4])*yd, col="gray40", xpd=TRUE)
            text(chr[[i]][3], u[3]-diff(u[3:4])*yd/2, uchr[i], col="white", xpd=TRUE)
            rect(chr[[i]][1]-.5, u[3], chr[[i]][2]+.5, u[4], xpd=TRUE)
        }

    } else {

        image.plot(y, pos, t(lod), yaxt="n", xlab=ylab, ylab="",
              zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n", axes=F)
        image(y, pos, t(lod), axes=FALSE, add=TRUE, pos, y, lod, xaxt="n", ylab=ylab, xlab="",
              zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n", axes=F)
        mtext(text=y.text, side=1, line=0.3, at=y, las=1, cex=0.8)
        title(ylab="Chromosome", mgp=c(2, 0, 0))
        u <- par("usr")
        xd <- 0.04
        for(i in seq(along=chr)) {
            rect(u[1] - diff(u[1:2])*xd, chr[[i]][1]-.5 , u[1], chr[[i]][2]+.5, col="gray40", xpd=TRUE)
            text(u[1]-diff(u[1:2])*xd/2, chr[[i]][3], uchr[i], col="white", xpd=TRUE)
            rect(u[1], chr[[i]][1]-.5 , u[2], chr[[i]][2]+.5, xpd=TRUE)
        }
    }
}

#-------------------------------------------------------------------#
# Example, i.e. your data goes here!
#-------------------------------------------------------------------#

# myqtls = runQTL(genotype = geno,
#     phenotype = metabolome_data_mixednorm,
#     marker_info = mrk,
#     permute = T,
#     pca = T,
#     permute_alpha = 0.1,
#     save_file = "./qtl_endometabolome_23042015/rqtls.rda")
