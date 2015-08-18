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
    # Run rQTL scanone analysis with optional methods,
    # permute_sig and pca
    genotype, # a n x m genotype matrix containing (eg 1,2) at 'n' genotype markers in 'm' strains
    phenotype, # a n x m phenotype matrix containing measurements of 'm' phenotypes in 'n' strains
    marker_info, # GRanges object with information about genetic markers (chromosome and location)
    permute = T, # compute significance of each QTL LOD by permutation
    pca = F, # maximize QTL detection by removing confounders using PCA
    permute_alpha = 0.05,
    save_file = ""){
  # subset genotype data on metabolite data. transpose it.
  if ( !sum(colnames(genotype)%in%rownames(phenotype))>0 ) {
    cat("ERROR:Strain names in genotype matrix (columns) do not match strain names in phenotype matrix (rows\n")
    return(NULL)
  }
  genotype_subset = t( genotype[ ,rownames( phenotype ) ] )
  # add chr num to markers
  genotype_subset = rbind( gsub('chr','', as.character( GenomicRanges::seqnames( marker_info[ colnames( genotype_subset ) ] ) ) ),
   genotype_subset )
  genotype_subset = cbind( c( NULL, rownames( genotype_subset ) ), genotype_subset )
  colnames( genotype_subset )[ 1 ] = 'id'
  # add id column to phen data
  phenotype = cbind( phenotype, rownames( phenotype )  )
  colnames( phenotype )[ dim( phenotype )[2] ] = 'id'
  # write rqtl files
  write.table( genotype_subset, file =  ".tmpgen", sep = ",", col.names = T, row.names = F, quote=F )
  write.table( phenotype, file =  ".tmpphen", sep = ",", col.names = T, row.names = F, quote=F )
  # read.cross to import data into rqtl
  cat("Creating rQTL cross object for genotype/phenotype data\n")
  genphen = try( qtl::read.cross( format = "csvs", ".", genfile = ".tmpgen" , phefile=  ".tmpphen", genotypes = c( "1","2" ) ) )
  if ( class(genphen)[1]!="try-error" ) {
    # clean up
    file.remove(c(".tmpgen",".tmpphen"))
  } else {
    cat("ERROR: Could not create rQTL cross object from genotype and phenotype data\n")
    return(NULL)
  }
  genphen = calc.genoprob( genphen,step = 0 )
  if (pca) {
    cat("Removing principal components to increase number of detected QTLs\n")
    pc_removed = 0
    # as a first approximation, use lod score > 2.5 as a "true" QTL
    # since running permutation is expensive
    phenotype = genphen$pheno[,colnames(genphen$pheno)!="id"]
    qtls = sum( unlist( mclapply( seq(1,dim(phenotype)[2]),function( i ) { sum(scanone( genphen, pheno.col = i )$lod >= 2.5) } ) ) )
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
    genphen_mod$pheno = phenotype_mod
    qtls_mod = sum( unlist( mclapply( seq(1,dim(phenotype)[2]),function( i ) { sum(scanone( genphen_mod, pheno.col = i )$lod >= 2.5) } ) ) )
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
      genphen_mod$pheno = phenotype_mod
      qtls_mod = sum( unlist( mclapply( seq(1,dim(phenotype)[2]),function( i ) { sum(scanone( genphen_mod, pheno.col = i )$lod >= 2.5) } ) ) )
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
    genphen$pheno = phenotype
    to_r$qtls = c( mclapply( seq(1,dim(phenotype)[2]),function( i ) { scanone( genphen, pheno.col = i ) } ) )
    names(to_r$qtls) = colnames(genphen$phen)
    to_r$pc_removed = pc_removed
    to_r$phenotype = phenotype
  } else {
    to_r = list()
    to_r$qtls = c( mclapply( seq(1,dim(genphen$pheno)[2]-1),function( i ) { scanone( genphen, pheno.col = i ) } ) )
    names(to_r$qtls) = colnames(genphen$pheno)[seq(1,dim(genphen$pheno)[2]-1)]
  }

  to_r$cross = genphen
  if (permute) {
    # permuations to determine qtl sig cutoff
    to_r$qtls_permuted = lapply(seq(1,dim(phenotype)[2]), function(i){scanone( genphen, pheno.col = i, n.perm = 1000, n.cluster = 20 )})
    names(to_r$qtls_permuted) = colnames(genphen$phen)
    to_r$sig_qtls = list()
    to_r$qtls_threshold = summary(to_r$qtls_permuted,permute_alpha)
    for ( i in names(to_r$qtls) ) {
      #phe_qtls = summary( to_r$qtls[[ i ]] , to_r$qtls_threshold[,i] )
      phe_qtls = summary(to_r$qtls[[ i ]], format="allpeaks", perms=to_r$qtls_permuted[[ i ]],
        alpha=permute_alpha, pvalues=TRUE)
      if ( dim( phe_qtls )[ 1 ] > 0 ) {
        # format table entry
        add_info = c()
        for ( n in rownames( phe_qtls ) ) {
          add_info = rbind( add_info, as.data.frame( ranges( marker_info[ n ] ) ) )
        }
        phe_qtls = cbind( i, phe_qtls[ add_info$names, ], add_info )
        phe_qtls$pos = NULL
        to_r$sig_qtls[[i]] =  phe_qtls
      }
    }
    to_r$permute_alpha = permute_alpha
  }
  if (save_file!="") {
    qtl_list = to_r
    save(qtl_list,file=save_file)
  }
  return(to_r)
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
