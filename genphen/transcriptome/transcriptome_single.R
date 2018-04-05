#! /usr/bin/env Rscript
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/transcriptome/transcriptome.R")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Must provide gene as first input arg and n resamples as second input arg.n", call.=FALSE)
}

gene = args[1]
resamples = as.integer(args[2])
print(paste("Gene: ", gene))
print(paste("Resamples: ", resamples))

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
# only keep the protein coding genes
tokeep = grep("^Y",colnames(pheno))
#pheno = pheno[,tokeep]
pheno = pheno[,grep(gene,colnames(pheno)]

qtl_f = file.path(BASEDIR,"qtl", paste("eQTL_",gene,".rda"))
if (!file.exists(qtl_f)) {
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
	pb = txtProgressBar(min = 0, max = dim(eQTL$qtls)[2], style = 3)
	eQTL$resamples = lapply(seq(1,dim(eQTL$qtls)[2]), function(i) {
		setTxtProgressBar(pb, i)
		try(qtl::scanone( cross, pheno.col = phes[i], n.perm = resamples, n.cluster = 24))
	})
	close(pb)
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
