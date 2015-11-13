#! /usr/bin/env Rscript
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/transcriptome/transcriptome.R")

#-------------------------------------------------------------------#
# GenPhen segregant Transcriptome Analysis
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "WTFPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"

# Import packages ---------------------------------------------------

# library( plotly )
# set_credentials_file( "scalefreegan", "5xyvvaleim" )
library( qtl )
library( GenomicRanges )
library( rtracklayer )
library( plyr )
library( parallel )
options("mc.cores"=24)
library( ggplot2 )
library( ggbio )
library( LSD )
library( DESeq2 )
library(dplyr)
library(reshape2)
library(pheatmap)

# source rQTL utilities
devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/QTL/rQTL.R")

# load Chenchen's processed transcriptome data
load("/g/steinmetz/project/GenPhen/data/3tagseq/counts.rda")

BASEDIR = "/g/steinmetz/brooks/genphen/transcriptome"
.plot = FALSE

#-------------------------------------------------------------------#
# Identify strains with low coverage
#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
# Correct
#-------------------------------------------------------------------#

# QC ---------------------------------------------------

if (.plot) {
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
	pdf(file=file.path(BASEDIR, "plots","transcriptome_biorep_reproducibility.pdf"))
		heatscatter(rep_mat_log[,1],rep_mat_log[,2],main="Transcript counts: Reproducibility across biological replicates", cor=T, xlab="Rep1, log2(transcript count)", ylab="RepX, log2(transcript count)")
	dev.off()

	# replicate perstrain
	pdf(file=file.path(BASEDIR, "plots","transcriptome_biorep_reproducibility_perstrain.pdf"))
	d = select(trx_df, strain, rep, value) %>%
		group_by(.,strain) %>%
	do({
		x = filter(.,rep=="biorep1")
		y = filter(.,rep=="biorep2")
		x.log2 = log2(x$value+1)
		y.log2 = log2(y$value+1)
		if (length(x.log2)==length(y.log2)) {
			data.frame(x.log2 = x.log2, y.log2 = y.log2)
		} else {
			data.frame()
		}
		})
	ggplot(d, aes(x.log2, y.log2)) +
		#geom_point() +
		stat_density2d(aes(fill = ..level..), geom="polygon") +
		facet_wrap(~ strain)
	dev.off()

	dcor = select(trx_df, strain, rep, value) %>%
		group_by(.,strain) %>%
		do({
				o1 = try(cor(log2(.$value[.$rep=="biorep1"]+1),log2(.$value[.$rep=="biorep2"]+1),use="pair"),silent=T)
				if (class(o1)=="try-error") {
					o1 = NA
				}
				return(data.frame(cor = o1))
			})
	 pdf(file=file.path(BASEDIR, "plots","transcriptome_biorep_reproducibility_perstrain_corhist.pdf"))
	 	hist(dcor$cor,100,xlim=c(0,1),main="Biorep correlation", xlab="Pearson correlation")
		hist(dcor$cor,100,xlim=c(0.9,1),main="Biorep correlation", xlab="Pearson correlation")
	 dev.off()




}

#-------------------------------------------------------------------#
# Rqtl
#-------------------------------------------------------------------#

strainRename = function(strains) {
	o = sapply(strains,function(strain){
		if (nchar(strain)==3) {
			strain = gsub("^X","0",strain)
		} else {
			strain = gsub("^X","",strain)
		}
		return(strain)
	})
	return(o)
}

f = file.path(BASEDIR, "data","trx_df.rda")
# This is a quick and dirty first shot at Rqtl
# Trx data needs more quality control before moving forward
if (!file.exists(f)) {
	# reformat data
	pb = txtProgressBar(min = 0, max = dim(normCounts)[1], style = 3)
	trx_df = lapply(seq(1,dim(normCounts)[1]), function(i) {
		setTxtProgressBar(pb, i)
		#print(i)
		d = normCounts[i,]
		n = do.call(rbind,strsplit(names(d),split="_"))
		colnames(n) = c("strain","rep")
		n[,"strain"] = strainRename(n[,"strain"])
		o = data.frame(
			name = tx_3utr[i,"Name"],
			type = tx_3utr[i,"type"],
			chr = tx_3utr[i,"seqnames"],
			start = tx_3utr[i,"start"],
			end = tx_3utr[i,"end"],
			strand = tx_3utr[i,"strand"],
			value = d
		)
		o = cbind(n,o)
		return(o)
	})
	close(pb)
	trx_df = do.call(rbind,trx_df)
	save(trx_df, file = f)
} else {
	load(f)
}

# load genotype and markers files
genotype_f = "/g/steinmetz/brooks/yeast/genomes/S288CxYJM789/genotypes_S288c_R64.rda"
load(genotype_f)

# reduce data set - take mean of bioreps
pheno = acast(trx_df, formula =strain~name, fun.aggregate = mean, value.var = "value")
# only take strains in genotype matrix
# take the shifted logarithm
pheno = log2(pheno[intersect(rownames(pheno),colnames(geno)),]+1)

qtl_f = file.path(BASEDIR,"qtl","eQTL.rda")
if (F) {
	cross_f = file.path(BASEDIR,"data","trx_cross.rda")
	if (!file.exists(cross_f)) {
		cross =	runQTL(
					genotype = geno,
					phenotype = pheno,
					marker_info = mrk,
					permute = F, # compute significance of each QTL LOD by permutation
			    pca = F, # maximize QTL detection by removing confounders using PCA
			    permute_alpha = 0.05,
			    save_file = qtl_f,
			    return_cross = TRUE, # just return cross object,
			    subset_genotype = F,
					estimate.map=FALSE,
					)
		save(cross, file = cross_f)
	} else {
		load(cross_f)
	}
	eQTL = list()
  # phenotype = cross$pheno[,colnames(cross$pheno)!="id",drop=F]
  # cross$pheno = phenotype
  # eQTL$qtls = c( lapply(seq(1,dim(phenotype)[2]),function( i ) { qtl::scanone( cross, pheno.col = i ) } ) )
	phes = which(colnames(cross$pheno)!="id")
	eQTL$qtls = qtl::scanone(cross, pheno.col = phes, method = "hk")
	eQTL$resamples = qtl::scanone( cross, pheno.col = phes, n.perm = 1000, n.cluster = 6)
	save(eQTL, file = qtl_f)
} else {
	load(qtl_f)
}

if (.plot) {
	# plot eqtls
	for (i in 3:dim(eQTL$qtls)[2]) {
		pdf(paste("/g/steinmetz/brooks/genphen/transcriptome/plots/lods/",colnames(eQTL$qtls)[i],".pdf",sep=""))
			try({
				tmp = tmp = eQTL$qtls[,1:3]
				tmp[,3] = eQTL$qtls[,i]
				plot(tmp, main = colnames(eQTL$qtls)[i], bandcol="gray90", ylab = "LOD")
			})
		dev.off()

		eQTL_cor = cor(eQTL$qtls[,3:dim(eQTL$qtls)[2]], use = "pair", method = "spearman")
		while(sum(is.na(eQTL_cor))>0) {
			tor = sort(apply(eQTL_cor,1,function(i){sum(is.na(i))}),decreasing=T)
			tor1 = which(rownames(eQTL_cor)==names(tor)[1])
			tor2 = which(colnames(eQTL_cor)==names(tor)[1])
			eQTL_cor = eQTL_cor[-tor1,-tor2]
		}
		pdf(paste("/g/steinmetz/brooks/genphen/transcriptome/plots/eQTL_cor.pdf"))
				pheatmap(eQTL_cor,breaks=seq(-1,1,length.out=100))
		dev.off()

		jpeg(paste("/g/steinmetz/brooks/genphen/transcriptome/plots/eQTL_cor.pdf"))
				tmp = tmp = eQTL$qtls[,1:3]
				tmp[,3] = eQTL$qtls[,i]
				plot(tmp, main = colnames(eQTL$qtls)[i], bandcol="gray90", ylab = "LOD")
		dev.off()
	}
}
