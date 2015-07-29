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

#-------------------------------------------------------------------#
# RF-QTL
#-------------------------------------------------------------------#

# Import packages ---------------------------------------------------

library( GenomicRanges )
library( rtracklayer )
library( plyr )
library( ColorPalettes )
#library( snow )
library( parallel )
options("mc.cores"=20)

# load RF functions ---------------------------------------------------

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
write.table( rfqtl_table, file =  "./qtl_endometabolome_23042015/rf_qtls_15062014.txt", sep = "\t", col.names = T, row.names = F, quote=F )


# save
save( rfqtl_table, endometabolome_rfqtls, file = "./qtl_endometabolome_23042015/rf_qtls.RData" )

# plot
# par(mfrow = c(2, 1))
# plot(sf["ARG",], type = "h", ylab = "adj. selection frequency", main = "RFSF, ARG")
# plot(sf["SUC",], type = "h", ylab = "adj. selection frequency", main = "RFSF, SUC")

# plot interaction differences across all tested interactions
# looking for very low/high differences
interaction_diffs = unlist(lapply(endometabolome_rfqtls,function(i){i$interaction[,"Difference"]}))
interaction_diffs_ecdf = ecdf( interaction_diffs )
sig_interactions = sort( interaction_diffs[ c( which( interaction_diffs < quantile(interaction_diffs_ecdf,.05) ), which( interaction_diffs > quantile(interaction_diffs_ecdf,.95) ) ) ] )
sig_interactions_separable = do.call( rbind,lapply( names( sig_interactions ), function(i) {
    info = mrk[strsplit( strsplit(i,"\\.")[[1]][2], ":")[[1]]]
    if ( length(runValue(seqnames(info))) > 1 ) {
        cbind(i,as.data.frame(info),sig_interactions[i])
      } else if ( ( start(ranges(info)[2]) - end(ranges(info)[1]) ) > 5000) {
        cbind(i,as.data.frame(info),sig_interactions[i])
        } else {
        return(NULL)
      }
  }) )
write.table( sig_interactions_separable, file =  "./qtl_endometabolome_23042015/sig_interactions_15062014.txt", sep = "\t", col.names = T, row.names = T, quote=F )
h = hist(interaction_diffs,100, plot=F)
cuts <- cut(h$breaks, c(-Inf,quantile(interaction_diffs_ecdf,.05),quantile(
  interaction_diffs_ecdf,.95),Inf))
pdf("./qtl_endometabolome_23042015/rf_qtl_interactions.pdf")
  plot(h, col=c("red","white","red")[cuts], main = "RF-QTL Pairwise Interactions, All metabolites", xlab = "Paired - Additive Importance")
dev.off()

pdf("./qtl_endometabolome_23042015/rf_qtl_correlation.pdf")
  tmp = do.call(rbind,lapply(endometabolome_rfqtls,function(i)i$sf))
  pheatmap(cor(t(tmp)),breaks=seq(-1,1,len=100))
dev.off()
