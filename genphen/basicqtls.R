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
library( randomForest )
#library( snow )
library( parallel )
options("mc.cores"=20) 

# Functions ---------------------------------------------------
source("/g/steinmetz/brooks/git/R-tools/quantile_normalize.R")

# Load data ---------------------------------------------------------
# NOTE: metabolome data was cleaned up. for strains < 10, a leading zero was added, e.g. 01B
metabolome_data = read.table( "./qtl_endometabolome_23042015/seg_endometabolome.txt", sep = "\t", header = T, row.names = 1, dec = "," )
load( "./qtl_endometabolome_23042015/geno_mrk.RData" )

# check distribution of the data

# apply shapiro-wilk test for normality

# raw
norm_test = p.adjust( apply( metabolome_data, 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 13/26 are not normal

# log10
norm_test_log = p.adjust( apply( log10( metabolome_data ), 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 8/26 are not normal

# zscore
norm_test_zscore = p.adjust( apply( apply( metabolome_data, 2, scale, center=TRUE, scale=TRUE ), 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 13/26 are not normal

# quantile normalize 
norm_test_quant = p.adjust( apply( quantile_normalize( metabolome_data ), 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 0/26 are not normal

# quantile normalize / zscore
norm_test_quant_zscore = p.adjust( apply( apply( quantile_normalize( metabolome_data ), 2, scale, center=TRUE, scale=TRUE ), 2, function(i) { shapiro.test( i )$p.value } ), "fdr" )
# 0/26 are not normal

# quantile normalize / zscore
# USE THIS!
metabolite_quant_zscore = apply( quantile_normalize( metabolome_data ), 2, scale, center=TRUE, scale=TRUE )
rownames( metabolite_quant_zscore ) = rownames( metabolome_data )

# as noted by Chenchen, log transform seems to take care of this


#-------------------------------------------------------------------#
# Rqtl 
#-------------------------------------------------------------------#

# the data format for rqtl really sucks...

# subset genotype data on metabolite data. transpose it.
geno_subset = t( geno[ ,rownames( metabolome_data ) ] )
# add chr num to markers
geno_subset = rbind( gsub('chr','', as.character( seqnames( mrk[ colnames( geno_subset ) ] ) ) ), geno_subset )
geno_subset = cbind( c( NULL, rownames( geno_subset ) ), geno_subset )
colnames( geno_subset )[ 1 ] = 'id'

# add id column to phen datayeas
# 4/5/2015 change to quant_zscore norm
# metabolome_data_w = cbind( metabolome_data, rownames( metabolome_data )  )
metabolome_data_w = cbind( metabolite_quant_zscore, rownames( metabolite_quant_zscore )  )
colnames( metabolome_data_w )[ dim( metabolome_data_w )[2] ] = 'id'

# write qtl files
write.table( geno_subset, file =  "./qtl_endometabolome_23042015/gen.csvs", sep = ",", col.names = T, row.names = F, quote=F )
write.table( metabolome_data_w, file =  "./qtl_endometabolome_23042015/phen.csvs", sep = ",", col.names = T, row.names = F, quote=F )

# read.cross to import data into rqtl
genphen = read.cross( format = "csvs", genfile = "./qtl_endometabolome_23042015/gen.csvs" , phefile=  "./qtl_endometabolome_23042015/phen.csvs", genotypes = c( "1","2" ) )

genphen = calc.genoprob( genphen,step = 0 )
qtls = c( mclapply( seq( 1,26 ),function( i ) { scanone( genphen, pheno.col = i ) } ) )
# permuations to determine qtl sig cutoff
qtls_permuted = c( lapply( seq( 1,26 ),function( i ) { scanone( genphen, pheno.col = i, n.perm = 1000, n.cluster = 20 ) } ) )

# Make table of putative QTLs for yeastmine
qtl_table = c()
permute_alpha = 0.05
for ( i in seq( 1, length(qtls) ) ) {
	phe = colnames( genphen$pheno )[ i ]
	qtl_threshold = summary( qtls_permuted[[ i ]], alpha = permute_alpha )[ 1 ]
	phe_qtls = summary( qtls[[ i ]] , qtl_threshold )
	if ( dim( phe_qtls )[ 1 ] > 0 ) {
		# format table entry
		add_info = c()
		for ( n in rownames( phe_qtls ) ) {
			add_info = rbind( add_info, as.data.frame( ranges( mrk[ n ] ) ) )
		}
		phe_qtls = cbind( phe, phe_qtls[ add_info$names, ], add_info )
		phe_qtls$pos = NULL
		qtl_table = rbind( qtl_table, phe_qtls )
	}
}
# sort by lod_score, phe, chr, start
qtl_table = arrange( qtl_table, desc( lod ), phe, chr, start )

# write table
write.table( qtl_table, file =  "./qtl_endometabolome_23042015/qtls_04052014.txt", sep = "\t", col.names = T, row.names = F, quote=F )


#-------------------------------------------------------------------#
# RF-QTL 
#-------------------------------------------------------------------#

# load RF functions
source( "./qtl_endometabolome_23042015/rfqtl_functions.R" )

#cl = makeSOCKcluster( 20 )
#clusterExport( cl, "rfsf" )

ff = function(index, data, x, ntree, corr, alpha ) {
	# actual data
  y = data[,index]
	rf = randomForest::randomForest(y = y, x = x, ntree = ntree)
	sf = rfsf( rf ) - corr
	# permuted background
	sf.null = numeric()
	for (i in 1:10) {
		rf.null = randomForest::randomForest( y = sample( y ), x = x, ntree = ntree )
		tmp = rfsf( rf.null ) - corr
		tmp[ tmp < 0 ] = 0
		sf.null = c( sf.null, tmp )
	}
	pnull = ecdf(sf.null)
	P = 1 - pnull(sf)
	Q = p.adjust(P, "fdr")
	to_r = list()
	to_r$qtls = sf[ which( Q <= alpha ) ]
	to_r$sf = sf
	to_r$Q = Q
	to_r$ntree = ntree
	to_r$alpha = alpha
	return( to_r )
}

# data
corr = estBias(metabolite_quant_zscore, 10000, verbose = F)
#endometabolome_rfqtls = parApply(cl, t( metabolome_data ), 1, ff, t( geno )[ rownames( metabolome_data),  ], 2000, corr, 0.05 )
endometabolome_rfqtls = mclapply(colnames( metabolite_quant_zscore ), ff, y = metabolite_quant_zscore[ ,i ], x = t( geno )[ rownames( metabolome_data),  ], ntree = 2000, corr = corr, alpha = 0.05 )

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

# normal
z = as.matrix( metabolome_data[ metabolome_data_clustx$order, metabolome_data_clusty$order ] )
rownames( z ) = NULL
colnames( z ) = NULL

data = list(
  list(
    z = z, 
    x = colnames( metabolome_data ), 
    y = rownames( metabolome_data ), 
    colorscale = "Jet",
    type = "heatmap"
  )
)

response = py$plotly(data, kwargs=list(filename="Endometabolome", fileopt="overwrite"))
url = response$url

# log
z_log = log( as.matrix( metabolome_data[ metabolome_data_clustx$order, metabolome_data_clusty$order ] ) )
rownames( z_log ) = NULL
colnames( z_log ) = NULL

data = list(
  list(
    z = z_log, 
    x = colnames( metabolome_data ), 
    y = rownames( metabolome_data ), 
    colorscale = "Jet",
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