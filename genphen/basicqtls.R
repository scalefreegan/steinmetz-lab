#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# Endometabolome QTL Analysis
#-------------------------------------------------------------------# 

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = ["Aaron Brooks"]
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"

# Import packages ---------------------------------------------------
library( plotly )
# set_credentials_file( "scalefreegan", "5xyvvaleim" )
library( qtl )
library( GenomicRanges )
library( rtracklayer )
library( plyr )
library( ColorPalettes )
#library( snow )
library( parallel )
options("mc.cores"=20)

library( ggplot2 ) 
library( ggbio )

# yeast genome sequence
library("BSgenome.Scerevisiae.UCSC.sacCer3")
genome = BSgenome.Scerevisiae.UCSC.sacCer3

# load GRanges file
load("/g/steinmetz/brooks/git/steinmetz-lab/general/s288c_GRanges.RData")

# Functions ---------------------------------------------------
source( "/g/steinmetz/brooks/git/R-tools/quantile_normalize.R" )
source( "/g/steinmetz/brooks/git/steinmetz-lab/genphen/metabnorm/normfunctions.R" )

# Load data ---------------------------------------------------------
# NOTE: metabolome data was cleaned up. for strains < 10, a leading zero was added, e.g. 01B
metabolome_data = read.table( "./qtl_endometabolome_23042015/seg_endometabolome.txt", sep = "\t", header = T, row.names = 1, dec = "," )
# log the data
metabolome_data = log10( metabolome_data )
metabolome_data_quant = t(quantile_normalize(t(metabolome_data)))

#              #
# USE THIS ONE
#              #
metabolome_data_mixednorm = t( normalizeMixed( t( metabolome_data ), cohort = NULL, batch = NULL ) )

load( "./qtl_endometabolome_23042015/geno_mrk.RData" )

# QC data ---------------------------------------------------

# check distribution of the data

# apply shapiro-wilk test for normality

# raw log10
norm_test = p.adjust( apply( metabolome_data, 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 8/26 are not normal

# zscore
norm_test_zscore = p.adjust( apply( apply( metabolome_data, 2, scale, center=TRUE, scale=TRUE ), 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 8/26 are not normal

# quantile normalize 
norm_test_quant = p.adjust( apply( metabolome_data_quant, 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 21/26 are not normal

# quantile normalize / zscore
norm_test_quant_zscore = p.adjust( apply( apply( metabolome_data_quant, 2, scale, center=TRUE, scale=TRUE ), 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 21/26 are not normal

# boxplot(t(metabolome_data),main="Strains, Log10(Value)",xlab="Strain",las=2, ylab = "Metabolite Level",ylim=c(0,5))
# boxplot(t(metabolome_data_quant),main="Strains, Quantile Normalized",xlab="Strain",las=2,ylab = "Metabolite Level", ylim=c(0,5))

# boxplot((metabolome_data),main="Metabolites, Log10(Value)",xlab="Metabolite",las=2, ylab = "Metabolite Level", ylim=c(0,5))
# boxplot((metabolome_data_quant),main="Metabolites, Quantile Normalized",xlab="Metabolite",ylab = "Metabolite Level", las=2,ylim=c(0,5))

# boxplot(t(metabolome_data_mixednorm),main="Strains, Mix-Norm",xlab="Strain",las=2,ylab = "Metabolite Level", ylim=c(0,5))
# boxplot((metabolome_data_mixednorm),main="Metabolites, Mix-Norm",xlab="Metabolite",las=2, ylab = "Metabolite Level", ylim=c(0,5))

# reduced metabolite correlations
# m_cor = cor(metabolome_data)[upper.tri(cor(metabolome_data))]
# m_cor_mixednorm = cor(metabolome_data_mixednorm)[upper.tri(cor(metabolome_data_mixednorm))]
#boxplot(list( m_cor,m_cor_mixednorm ),main="Metabolite Correlation",xlab="Normalization method",las=1, ylab = "Metabolite Correlation, Pearson", ylim=c(-1,1),names=list("Log10","Mixed Model"))


# QC genotype ---------------------------------------------------

# asses allele frequencies

g_1 = apply( geno, 1, function(x) { sum(x==1) })
g_2 = apply( geno, 1, function(x) { sum(x==2) })

# ggbio way
p <- ggbio( trackWidth = 1, buffer = 0, radius = 20 ) + 
  circle(yeast_gr, geom = "ideo", fill = "gray70") +
  circle(yeast_gr, geom = "scale", size = 2) +
  circle(yeast_gr, geom = "text", aes(label = seqnames), vjust = 0, size = 3)

# circlize way
df = as.data.frame(yeast_gr)
mrk_df = as.data.frame(mrk)
chr_name_ttable = levels(df$seqnames)
names(chr_name_ttable) = 

circos.par(gap.degree = 2)
circos.genomicInitialize(df, track.height = .1)
circos.genomicTrackPlotRegion(ylim = c(0, 1),
bg.col = rep("#838B8B", dim(df)[1]),
bg.border = NA, track.height = 0.05)

#-------------------------------------------------------------------#
# Rqtl 
#-------------------------------------------------------------------#

genotype = geno
phenotype = metabolome_data_mixednorm
marker_info = mrk

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
  genotype_subset = t( genotype[ ,rownames( metabolome_data ) ] )
  # add chr num to markers
  genotype_subset = rbind( gsub('chr','', as.character( seqnames( marker_info[ colnames( genotype_subset ) ] ) ) ), genotype_subset )
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
  genphen = try( read.cross( format = "csvs", ".", genfile = ".tmpgen" , phefile=  ".tmpphen", genotypes = c( "1","2" ) ) )
  if ( class(genphen)!="try-error" ) {
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
    qtls = sum( unlist( mclapply( seq( 1,26 ),function( i ) { sum(scanone( genphen, pheno.col = i )$lod >= 2.5) } ) ) )
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
    qtls_mod = sum( unlist( mclapply( seq( 1,26 ),function( i ) { sum(scanone( genphen_mod, pheno.col = i )$lod >= 2.5) } ) ) )
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
      qtls_mod = sum( unlist( mclapply( seq( 1,26 ),function( i ) { sum(scanone( genphen_mod, pheno.col = i )$lod >= 2.5) } ) ) )
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
    to_r$qtls = c( mclapply( seq( 1,26 ),function( i ) { scanone( genphen, pheno.col = i ) } ) )
    names(to_r$qtls) = colnames(genphen$phen)
    to_r$pc_removed = pc_removed
    to_r$phenotype = phenotype
  } else {
    to_r$qtls = c( mclapply( seq( 1,26 ),function( i ) { scanone( genphen, pheno.col = i ) } ) )
    names(to_r$qtls) = colnames(genphen$phen)
  }

  to_r$cross = genphen
  if (permute) {
    # permuations to determine qtl sig cutoff
    to_r$qtls_permuted = c( lapply( seq( 1,26 ),function( i ) { scanone( genphen, pheno.col = i, n.perm = 1000, n.cluster = 20 ) } ) )
    names(to_r$qtls_permuted) = colnames(genphen$phen)
    to_r$sig_qtls = list()
    to_r$qtls_threshold = list()
    for ( i in seq( 1, length(to_r$qtls) ) ) {
      phe = colnames( genphen$pheno )[ i ]
      to_r$qtls_threshold[[phe]] = summary( qtls_permuted[[ i ]], alpha = permute_alpha )[ 1 ]
      phe_qtls = summary( to_r$qtls[[ i ]] , to_r$qtl_threshold[[phe]] )
      if ( dim( phe_qtls )[ 1 ] > 0 ) {
        # format table entry
        add_info = c()
        for ( n in rownames( phe_qtls ) ) {
          add_info = rbind( add_info, as.data.frame( ranges( marker_info[ n ] ) ) )
        }
        phe_qtls = cbind( phe, phe_qtls[ add_info$names, ], add_info )
        phe_qtls$pos = NULL
        to_r$sig_qtls[[phe]] =  phe_qtls
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

myqtls = runQTL(genotype = geno, 
    phenotype = metabolome_data_mixednorm,
    marker_info = mrk,
    permute = T, 
    pca = T, 
    permute_alpha = 0.1,
    save_file = "./qtl_endometabolome_23042015/rqtls.rda")

# # Make table of putative QTLs for yeastmine
# qtl_table = c()
# permute_alpha = 0.05
# for ( i in seq( 1, length(qtls) ) ) {
# 	phe = colnames( genphen$pheno )[ i ]
# 	qtl_threshold = summary( qtls_permuted[[ i ]], alpha = permute_alpha )[ 1 ]
# 	phe_qtls = summary( qtls[[ i ]] , qtl_threshold )
# 	if ( dim( phe_qtls )[ 1 ] > 0 ) {
# 		# format table entry
# 		add_info = c()
# 		for ( n in rownames( phe_qtls ) ) {
# 			add_info = rbind( add_info, as.data.frame( ranges( mrk[ n ] ) ) )
# 		}
# 		phe_qtls = cbind( phe, phe_qtls[ add_info$names, ], add_info )
# 		phe_qtls$pos = NULL
# 		qtl_table = rbind( qtl_table, phe_qtls )
# 	}
# }
# # sort by lod_score, phe, chr, start
# qtl_table = arrange( qtl_table, desc( lod ), phe, chr, start )

# # write table
# write.table( qtl_table, file =  "./qtl_endometabolome_23042015/qtls_04052014.txt", sep = "\t", col.names = T, row.names = F, quote=F )


#-------------------------------------------------------------------#
# RF-QTL 
#-------------------------------------------------------------------#

# load RF functions

# DO NOT USE THIS AS OF 5/6/2015. Massive memory leak
#library( parallelRandomForest )
# randomForestSRC package contains fcns for interactions
library(randomForestSRC)

#source( "./qtl_endometabolome_23042015/rfqtl_functions.R" )
# load functions rewritten for rfsrc 
source( "./qtl_endometabolome_23042015/rfsrcqtl_functions.R" )

#cl = makeSOCKcluster( 20 )
#clusterExport( cl, "rfsf" )


ff = function(index, y, x, ntree, corr, alpha,interaction=F,... ) {
	# one metabolite data
  y = y[,index]
  # composite data/genotype data frame
  data = as.data.frame( cbind( y[rownames(x)] , x ) )
  rf = randomForestSRC::rfsrc( V1 ~ ., data=data, ntree=ntree, var.used="all.trees",forest=TRUE)
	sf = rfsf( rf ) 
  sf = sf - corr[ names(sf) ]
	# permuted background
	sf.null = numeric()
	for (i in 1:10) {
    permuted_data = data
    permuted_data[,1] = sample( permuted_data[,1] )
		rf.null = randomForestSRC::rfsrc( V1 ~ ., data=permuted_data, ntree=ntree, var.used="all.trees")
		tmp = rfsf( rf.null ) 
    tmp = tmp - corr[ names(tmp) ]
		tmp[ tmp < 0 ] = 0
		sf.null = c( sf.null, tmp )
	}
	pnull = ecdf(sf.null)
	P = 1 - pnull(sf)
	Q = p.adjust(P, "fdr")
  names(Q) = names(sf)
	to_r = list()
  to_r$rf = rf
	to_r$qtls = sf[ which( Q <= alpha ) ]
	to_r$sf = sf
	to_r$Q = Q
	to_r$ntree = ntree
	to_r$alpha = alpha
  if (interaction) {
    try({
    # compute pairwise interaction
    # use method="maxsubtree"
    # to reduce computation time only consider only QTLs or top 5 predictors
    if ( length( to_r$qtls )>0 ){
      to_test = to_r$qtls
    } else {
      to_test = names( sort( to_r$Q ) )[1:5] 
    }
    to_r$interaction = randomForestSRC::find.interaction( rf,  xvar.names = to_test, method="vimp" )
    })
  }
	return( to_r )
}

# data
# corr = estBias( metabolome_data_mixednorm, 10000, verbose = F )
x = t( geno )[ rownames( metabolome_data_mixednorm ),  ]
# do 10000
if (file.exists("./qtl_endometabolome_23042015/rf_mrkbias.RData") ) {
  load("./qtl_endometabolome_23042015/rf_mrkbias.RData")
} else {
  corr = estBias( x, 10000, verbose = T )
  save(corr,file="./qtl_endometabolome_23042015/rf_mrkbias.RData")
}
#endometabolome_rfqtls = parApply(cl, t( metabolome_data ), 1, ff, t( geno )[ rownames( metabolome_data),  ], 2000, corr, 0.05 )
endometabolome_rfqtls = lapply( colnames( metabolome_data_mixednorm ), ff, y = metabolome_data_mixednorm, x = x, ntree = 2000, corr = corr, alpha = 0.05, interaction = T )
names( endometabolome_rfqtls ) = colnames( metabolome_data_mixednorm )

rfqtl_table = c()
for ( i in names(endometabolome_rfqtls ) ) {
	tmp_qtls = endometabolome_rfqtls[[i]]$qtls
	if ( length( tmp_qtls ) > 0 ) {
		# format table entry
		add_info = c()
		for ( n in names( tmp_qtls ) ) {
			chr = gsub( "chr", "", as.character( seqnames( mrk[ n ] ) ) )
			add_info = rbind( add_info, cbind( chr, as.data.frame( ranges( mrk[ n ] ) ) ) )
		}
		phe = i
		selection_freq = tmp_qtls[ add_info$names ]
		tmp_qtls = cbind( phe, selection_freq, add_info )
		rfqtl_table = rbind( rfqtl_table, tmp_qtls )
	}
}
# sort by lod_score, phe, chr, start
rfqtl_table = arrange( rfqtl_table, phe, desc( selection_freq ), chr, start )

# write table
#write.table( rfqtl_table, file =  "./qtl_endometabolome_23042015/rf_qtls.txt", sep = "\t", col.names = T, row.names = F, quote=F )
write.table( rfqtl_table, file =  "./qtl_endometabolome_23042015/rf_qtls_04052014.txt", sep = "\t", col.names = T, row.names = F, quote=F )


# save
save( rfqtl_table, endometabolome_rfqtls, file = "./qtl_endometabolome_23042015/rf_qtls.RData" )

# plot
# par(mfrow = c(2, 1))
# plot(sf["ARG",], type = "h", ylab = "adj. selection frequency", main = "RFSF, ARG")
# plot(sf["SUC",], type = "h", ylab = "adj. selection frequency", main = "RFSF, SUC")


#-------------------------------------------------------------------#
# Plot
#-------------------------------------------------------------------#

# Cluster metabolome data ---------------------------------------------------------
metabolome_data_clustx = hclust( dist( metabolome_data ) )
metabolome_data_clusty = hclust( dist( t( metabolome_data ) ) )

# Plot metabolite data

py = plotly()

#
# raw data
#

# log
y = metabolome_data_clusty$labels[metabolome_data_clusty$order]
x = metabolome_data_clustx$labels[metabolome_data_clustx$order]
z = as.matrix( metabolome_data[ x, y ] )
rownames( z ) = NULL
colnames( z ) = NULL

data = list(
  list(
    z = z, 
    x = y, 
    y = x, 
    colorscale = "Jet",
    zauto = F,
    zmin = 0,
    zmax = 5, 
    type = "heatmap"
  )
)

response = py$plotly(data, kwargs=list(filename="Endometabolome_log", fileopt="overwrite"))
url = response$url


data = lapply( colnames( metabolite_quant_zscore ), function(i) {
    out = list()
    out$x = metabolite_quant_zscore[ ,i ]
    out$opacity = 0.75
    out$type = "histogram"
    out$name = i
    return(out)
    } ) 
layout <- list( 
  barmode = "overlay",
  title = "Distribution of Endometabolite Measurements (Quantile normalized, Z-score)", 
  xaxis = list(title = "Metabolite Level (Z-score)"), 
  yaxis = list(title = "Count") 
  )
response <- py$plotly(data, kwargs=list(layout=layout, filename="endometabolome histogram", fileopt="overwrite"))
url <- response$url

# mixed model normalized
y = metabolome_data_clusty$labels[metabolome_data_clusty$order]
x = metabolome_data_clustx$labels[metabolome_data_clustx$order]
z = as.matrix( metabolome_data_mixednorm[ x, y ] )
rownames( z ) = NULL
colnames( z ) = NULL

data = list(
  list(
    z = z, 
    x = y, 
    y = x, 
    colorscale = "Jet",
    zauto = F,
    zmin = 0,
    zmax = 5, 
    type = "heatmap"
  )
)

response = py$plotly(data, kwargs=list(filename="Endometabolome_mixednorm", fileopt="overwrite"))
url = response$url


data = lapply( colnames( metabolite_quant_zscore ), function(i) {
    out = list()
    out$x = metabolite_quant_zscore[ ,i ]
    out$opacity = 0.75
    out$type = "histogram"
    out$name = i
    return(out)
    } ) 
layout <- list( 
  barmode = "overlay",
  title = "Distribution of Endometabolite Measurements (Quantile normalized, Z-score)", 
  xaxis = list(title = "Metabolite Level (Z-score)"), 
  yaxis = list(title = "Count") 
  )
response <- py$plotly(data, kwargs=list(layout=layout, filename="endometabolome histogram", fileopt="overwrite"))
url <- response$url

#
# correlations
#

# metabolites
metabolite_cor = cor( metabolome_data, method="spearman" )
metabolite_cor_clustx = hclust( dist( metabolite_cor ) )
metabolite_cor_clusty = hclust( dist( t( metabolite_cor ) ) )

z = as.matrix( metabolite_cor[ metabolite_cor_clustx$order, metabolite_cor_clusty$order ] )
rownames( z ) = NULL
colnames( z ) = NULL

data = list(
  list(
    z = z, 
    x = colnames( metabolite_cor ), 
    y = rownames( metabolite_cor ), 
    colorscale = "Jet",
    zauto = F,
    zmin = -1,
    zmax = 1, 
    type = "heatmap"
  )
)

response = py$plotly(data, kwargs=list(filename="Endometabolome, Spearman Correlation", fileopt="overwrite"))
url = response$url

# strains
strain_cor = cor( t( metabolome_data ), method="spearman" )
strain_cor_clustx = hclust( dist( strain_cor ) )
strain_cor_clusty = hclust( dist( t( strain_cor ) ) )

z = as.matrix( strain_cor[ strain_cor_clustx$order, strain_cor_clusty$order ] )
rownames( z ) = NULL
colnames( z ) = NULL

data = list(
  list(
    z = z, 
    x = colnames( strain_cor ), 
    y = rownames( strain_cor ), 
    colorscale = "Jet",
    zauto = F,
    zmin = -1,
    zmax = 1, 
    type = "heatmap"
  )
)

response = py$plotly(data, kwargs=list(filename="Strains, Spearman Correlation", fileopt="overwrite"))
url = response$url

#-------------------------------------------------------------------#
# Analyse QTL effects 
#-------------------------------------------------------------------#

# endometabolome rQTLs
#   phe chr      lod  start    end width     names
# 1 ARG  13 4.705359  44287  44503   217  mrk_9421
# 2 SUC  16 4.210699 922017 922017     1 mrk_13643
# 3 LYS  15 4.136895  70894  70894     1 mrk_11483
# 4 GLN  15 3.752371 485898 485898     1 mrk_12115

# endometabolome rf_QTLs
#    phe selection_freq chr  start    end width     names
# 1  ARG    0.008146087  13  44287  44503   217  mrk_9421
# 2  ARG    0.006678967  13  43096  43166    71  mrk_9420
# 3  ARG    0.005760963  16 699992 700480   489 mrk_13402
# 4  ARG    0.004877374  13  63358  63947   590  mrk_9437
# 5  ARG    0.004843013  13  56748  56748     1  mrk_9435
# 6  ASN    0.008546134  02 295710 296082   373   mrk_581
# 7  GLU    0.006683078  15 493331 493377    47 mrk_12139
# 8  GLN    0.006000438  15 485898 485898     1 mrk_12115
# 9  MET    0.006332202  15 493331 493377    47 mrk_12139
# 10 TYR    0.005174024  15 480144 480336   193 mrk_12104

#
# 1. ARG
#
mrks = rfqtl_table[ rfqtl_table$phe == "ARG", "names" ]
genotype = t( geno[ mrks, rownames( metabolome_data ) ] )
arg_data = cbind( metabolome_data[ , "ARG", drop = F ], genotype )

x_1 = c( )
x_2 = c( )
y_1 = c( )
y_2 = c( )
for (i in mrks ) {
	x_1 = c( x_1, rep( i, sum( arg_data[ , i ] == 1 ) ) )
	x_2 = c( x_2, rep( i, sum( arg_data[ , i ] == 2 ) ) )
	y_1 = c( y_1, arg_data[ , "ARG" ][arg_data[ , i ] == 1] )
  y_2 = c( y_2, arg_data[ , "ARG" ][arg_data[ , i ] == 2] )
}

py = plotly()
data = list(
  list(
    y = y_1, 
    x = x_1,
    boxpoints = "all", 
    jitter = 0.3, 
    pointpos = 0, 
    type = "box",
    name = "S288c"
  ),
  list(
    y = y_2, 
    x = x_2,
    boxpoints = "all", 
    jitter = 0.3, 
    pointpos = 0, 
    type = "box",
    name = "YJM789"
  )
) 
layout = list(
  title = "ARG-related QTLs", 
  boxmode="group",
  xaxis = list(
    title ="QTL"
    ),
  yaxis = list(
    title = "ARG Level (unknown units)"
    )
  )
response = py$plotly(data, kwargs=list(layout=layout, filename="ARG QTL, Arg level", fileopt="overwrite"))
url = response$url

#
# 2. SUC
#
genotype = geno[ "mrk_13643", rownames( arg_data ) ]
arg_data = cbind( metabolome_data[ , "SUC", drop = F ], genotype )

py = plotly()
data = list(
  list(
    y = arg_data[ arg_data[ , "genotype" ] == 1, "SUC" ], 
    boxpoints = "all", 
    jitter = 0.3, 
    pointpos = 0, 
    type = "box",
    name = "S288c"
  ),
  list(
    y = arg_data[ arg_data[ , "genotype" ] == 2, "SUC" ], 
    boxpoints = "all", 
    jitter = 0.3, 
    pointpos = 0, 
    type = "box",
    name = "YJM789"
  )
) 
layout = list(title = "Suc level by variant at mrk_13643")
response = py$plotly(data, kwargs=list(layout=layout, filename="SUC QTL, Suc level", fileopt="overwrite"))
url = response$url

#-------------------------------------------------------------------#
# Load genomes 
#-----------------------------------------------------------------