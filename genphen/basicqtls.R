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
library( snow )

# Load data ---------------------------------------------------------
# NOTE: metabolome data was cleaned up. for strains < 10, a leading zero was added, e.g. 01B
metabolome_data = read.table( "./qtl_endometabolome_23042015/seg_endometabolome.txt", sep = "\t", header = T, row.names = 1, dec = "," )
load( "./qtl_endometabolome_23042015/geno_mrk.RData" )

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
metabolome_data_w = cbind( metabolome_data, rownames( metabolome_data )  )
colnames( metabolome_data_w )[ dim( metabolome_data_w )[2] ] = 'id'

# write qtl files
write.table( geno_subset, file =  "./qtl_endometabolome_23042015/gen.csvs", sep = ",", col.names = T, row.names = F, quote=F )
write.table( metabolome_data_w, file =  "./qtl_endometabolome_23042015/phen.csvs", sep = ",", col.names = T, row.names = F, quote=F )

# read.cross to import data into rqtl
genphen = read.cross( format = "csvs", genfile = "./qtl_endometabolome_23042015/gen.csvs" , phefile=  "./qtl_endometabolome_23042015/phen.csvs", genotypes = c( "1","2" ) )

genphen = calc.genoprob( genphen,step = 0 )
qtls = c( lapply( seq( 1,26 ),function( i ) { scanone( genphen, pheno.col = i ) } ) )
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
write.table( qtl_table, file =  "./qtl_endometabolome_23042015/qtls.txt", sep = "\t", col.names = T, row.names = F, quote=F )

#-------------------------------------------------------------------#
# RF-QTL 
#-------------------------------------------------------------------#

# load RF functions
source( "./qtl_endometabolome_23042015/rfqtl_functions.R" )

cl = makeSOCKcluster( 20 )
clusterExport( cl, "rfsf" )

ff = function(y, x, ntree ) {
	rf = randomForest::randomForest(y = y, x = x, ntree = ntree)
	sf = rfsf( rf )
}

# data
x = t( geno )[ rownames( metabolome_data),  ]
y = t( metabolome_data )
corr = estBias(x, 10000, verbose = F)
sf = parApply(cl, y, 1, ff, x, 1000 )
sf = t(sf - corr)

# single core way
y = metabolome_data[,"ARG"]
rf = randomForest(y = y, x = x, ntree = 2000)
sf = rfsf(rf)
sf.corr = sf - corr
sf.corr[sf.corr < 0] = 0

# sig
sf.null = numeric()
for (i in 1:10) {
	rf.null = randomForest(y = sample(y), x = x, ntree = 2000)
	tmp = rfsf(rf.null) - corr
	tmp[tmp < 0] = 0
	sf.null = c(sf.null, tmp)
	print(i)
}
pnull = ecdf(sf.null)
P = 1 - pnull(sf)
Q = p.adjust(P, "fdr")
names( which( sf>=Q ) )


# plot
par(mfrow = c(2, 1))
plot(sf["ARG",], type = "h", ylab = "adj. selection frequency", main = "RFSF, ARG")
plot(sf["SUC",], type = "h", ylab = "adj. selection frequency", main = "RFSF, SUC")


#-------------------------------------------------------------------#
# Plot
#-------------------------------------------------------------------#

# Cluster metabolome data ---------------------------------------------------------
metabolome_data_clustx = hclust( dist( metabolome_data ) )
metabolome_data_clusty = hclust( dist( t( metabolome_data ) ) )

# Plot metabolite data

py <- plotly()

#
# raw data
#

# normal
z = as.matrix( metabolome_data[ metabolome_data_clustx$order, metabolome_data_clusty$order ] )
rownames( z ) = NULL
colnames( z ) = NULL

data <- list(
  list(
    z = z, 
    x = colnames( metabolome_data ), 
    y = rownames( metabolome_data ), 
    colorscale = "Jet",
    type = "heatmap"
  )
)

response <- py$plotly(data, kwargs=list(filename="Endometabolome", fileopt="overwrite"))
url <- response$url

# log
z_log = log( as.matrix( metabolome_data[ metabolome_data_clustx$order, metabolome_data_clusty$order ] ) )
rownames( z_log ) = NULL
colnames( z_log ) = NULL

data <- list(
  list(
    z = z_log, 
    x = colnames( metabolome_data ), 
    y = rownames( metabolome_data ), 
    colorscale = "Jet",
    type = "heatmap"
  )
)

response <- py$plotly(data, kwargs=list(filename="Endometabolome_log", fileopt="overwrite"))
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

data <- list(
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

response <- py$plotly(data, kwargs=list(filename="Endometabolome, Spearman Correlation", fileopt="overwrite"))
url <- response$url

# strains
strain_cor = cor( t( metabolome_data ), method="spearman" )
strain_cor_clustx = hclust( dist( strain_cor ) )
strain_cor_clusty = hclust( dist( t( strain_cor ) ) )

z = as.matrix( strain_cor[ strain_cor_clustx$order, strain_cor_clusty$order ] )
rownames( z ) = NULL
colnames( z ) = NULL

data <- list(
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

response <- py$plotly(data, kwargs=list(filename="Strains, Spearman Correlation", fileopt="overwrite"))
url <- response$url

#-------------------------------------------------------------------#
# Analyse QTL effects 
#-------------------------------------------------------------------#

# endometabolome QTLs
#   phe chr      lod  start    end width     names
# 1 ARG  13 4.705359  44287  44503   217  mrk_9421
# 2 SUC  16 4.210699 922017 922017     1 mrk_13643
# 3 LYS  15 4.136895  70894  70894     1 mrk_11483
# 4 GLN  15 3.752371 485898 485898     1 mrk_12115

#
# 1. ARG
#
genotype = t( geno[ c("mrk_9420","mrk_9421","mrk_9437","mrk_13402"), rownames( metabolome_data ) ] )
genotype[ genotype == 1 ] = 
arg_data = cbind( metabolome_data[ , "ARG", drop = F ], genotype )

x_1 = c( )
x_2 = c( )
y_1 = c( )
y_2 = c( )
for (i in c("mrk_9420","mrk_9421","mrk_9437","mrk_13402") ) {
	x_1 = c( x_1, rep( i, sum( arg_data[ , i ] == 1 ) ) )
	x_2 = c( x_2, rep( i, sum( arg_data[ , i ] == 2 ) ) )
}

py <- plotly()
data <- list(
  list(
    y = arg_data[ arg_data[ , "genotype" ] == 1, "ARG" ], 
    x = x
    boxpoints = "all", 
    jitter = 0.3, 
    pointpos = 0, 
    type = "box",
    name = "S288c"
  ),
  list(
    y = arg_data[ arg_data[ , "genotype" ] == 2, "ARG" ], 
    x = x
    boxpoints = "all", 
    jitter = 0.3, 
    pointpos = 0, 
    type = "box",
    name = "YJM789"
  )
) 
layout <- list(title = "Arg level by variant at mrk_9421")
response <- py$plotly(data, kwargs=list(layout=layout, filename="ARG QTL, Arg level", fileopt="overwrite"))
url <- response$url

#
# 2. SUC
#
genotype = geno[ "mrk_13643", rownames( arg_data ) ]
arg_data = cbind( metabolome_data[ , "SUC", drop = F ], genotype )

py <- plotly()
data <- list(
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
layout <- list(title = "Suc level by variant at mrk_13643")
response <- py$plotly(data, kwargs=list(layout=layout, filename="SUC QTL, Suc level", fileopt="overwrite"))
url <- response$url

#-------------------------------------------------------------------#
# Load genomes 
#-----------------------------------------------------------------