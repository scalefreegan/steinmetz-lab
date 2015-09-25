#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# GenPhen segregant Transcriptome Analysis
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
library( parallel )
options("mc.cores"=20)
library( ggplot2 )
library( ggbio )
library( LSD )

# load Chenchen's processed transcriptome data
load("/g/steinmetz/project/GenPhen/data/3tagseq/counts.rda")

#-------------------------------------------------------------------#
# Identify strains with low coverage
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
# Correct
#-------------------------------------------------------------------#

# QC ---------------------------------------------------

#-------------------------------------------------------------------#
# Rqtl
#-------------------------------------------------------------------#

# check quant across bioreps
ids =  sapply(colnames(normCounts),function(i)strsplit(i,"_")[[1]][1])
ids_table = table( ids )

rep_mat = do.call( rbind, lapply( names( ids_table ), function( i ) {
	ind = which( ids == i )
	if ( length( ind ) > 2 ) {
		return( do.call( rbind, lapply( seq( 2:length( ind ) ), function(j) { normCounts[ ,c( ind[1], ind[j] ) ] } ) ) )
	} else if (length( ind ) < 2) {
		return( NULL )
	} else{
		return( normCounts[ , ind ] )
	}
} ) )
colnames(rep_mat) = c("Biorep_1","Biorep_X")
rep_mat_log = log2(rep_mat[rowSums(rep_mat>=1)>1,])
pdf(file="./qtl_transcriptome/transcriptcount_bioreps.pdf")
	heatscatter(rep_mat_log[,1],rep_mat_log[,2],main="Transcript counts: Reproducibility across biological replicates", cor=T, xlab="Rep1, log2(transcript count)", ylab="RepX, log2(transcript count)")
dev.off()
