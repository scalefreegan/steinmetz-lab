# ! /usr/bin/env Rscript
# designed to be saved as .Rprofile in working dir
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/growth/growthQTL.R")
#-------------------------------------------------------------------#
# Detect growthQTLs for GenPhen Data from
# growth data from Nicole (metabolomics experiments)
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "GPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.plot = FALSE

# Import packages ---------------------------------------------------
library(ggplot2)
#library(plotly)
library(plyr)
library(dplyr)
library(reshape2)
library(LSD)
library(qtl)
library(pheatmap)
library(parallel)
options(mc.cores = 24)
library(snow)
#library(clustQTL)
#devtools::install_github("scalefreegan/steinmetz-lab/clustQTL")
library(reshape2)
library(funqtl)
library(stringr)


strainRename = function(strains) {
	o = sapply(strains,function(strain){
		if (nchar(strain)==3) {
			strain = gsub("^X","0",strain)
		} else if (nchar(strain)==2) {
			strain = paste("0", strain, sep="")
		} else {
			strain = gsub("^X","",strain)
		}
		if (length(grep("\\.",strain))>0) {
			rep = as.numeric(strsplit(strsplit(strain, split = "_")[[1]][2], split = "\\.")[[1]][1]) + 1
			strain = paste(strsplit(strain, split = "_")[[1]][1], rep, sep = "_")
		}
		return(strain)
	})
	return(o)
}

# source rQTL utilities
devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/QTL/rQTL.R")

endo_f = "/g/steinmetz/project/GenPhen/data/endometabolome/data/endometabolite_full_12102015.rda"

# load complete endometabolite file
load(endo_f)

# load genotype and markers files
genotype_f = "/g/steinmetz/brooks/yeast/genomes/S288CxYJM789/genotypes_S288c_R64.rda"
load(genotype_f)

#-------------------------------------------------------------------#
# cellconc with funQTL - MLOD/SLOD
# No smoothing / PCA
#
#-------------------------------------------------------------------#

cc = filter(endometabolite, metabolite == "AKG", time_format == "absolute" ) %>%
     select(., strain, replicate, time, cellconc) %>%
     group_by(., strain, time) %>%
     summarise(., mean(cellconc))

cc = dcast(cc, strain~time, value.var = "mean(cellconc)")
rownames(cc) = sapply(as.character(cc$strain),strainRename)
cc = cc[,2:dim(cc)[2]]
# LOG2 transform to make more normal
cc = log2(cc)


if (FALSE) {
#if (!file.exists(f)) {
	# group phenotype timepoints

	cross_f = "/g/steinmetz/brooks/genphen/growth/cross_template.rda"
	#cross_f_abs = "/g/steinmetz/brooks/genphen/metabolome/cross_template_abs.rda"
	if (!file.exists(cross_f)) {
		cross =	runQTL(
					genotype = geno,
					phenotype = gr,
					marker_info = mrk,
					return_cross = TRUE,
					estimate.map = FALSE
					)
		save(cross, file = cross_f)
	} else {
		load(cross_f)
	}

	cross = calc.genoprob(cross, step = 0)
	# last phenotype column is the "id" tag
	pcols = seq(1, dim(cross$pheno)[2] - 1)
 	qtls = scanone(cross, method = "hk", pheno.col = pcols)
	eff = geteffects(cross, pheno.col = pcols)
	qtls_alt = scanoneF(cross, pheno.cols = pcols, method = "hk")
	# calc permutation threshold
	permout = scanoneF(cross, pheno.cols = pcols,
	                    method = "hk", n.perm = 1000, n.cluster = 20, verbose = F)
	# identify chrs with slod/lod above permute val
	chrs = unique(c(
		qtls_alt[qtls_alt[,"slod"]>=summary(permout)["5%","slod"],"chr"],
		qtls_alt[qtls_alt[,"mlod"]>=summary(permout)["5%","mlod"],"chr"]
		)
	)
	qtl_intervals = list()
	if (length(chrs)>0) {
		for (i in chrs) {
			qtl_intervals[[as.character(i)]] = mrk[rownames(bayesint(qtls_alt, chr = str_pad(i, 2, pad = "0"), prob=0.95, lodcolumn=1))]
		}
	}
	ccQTL = list(qtls = qtls, eff = eff,
				qtls_alt = qtls_alt, qtl_intervals = qtl_intervals,
				permout = permout, pcols = pcols)
	save(ccQTL, file = f)
	# plot stuff

	# effect plots
	pdf("/g/steinmetz/brooks/genphen/growth/plots/ccQTL_effects.pdf")
		plotlod(ccQTL$qtls, ccQTL$eff, ccQTL$pcols, gap=25, ylab="Time")
		mtext("ccQTL", side = 3)
	dev.off()

	m_slod = max(ccQTL$qtls_alt[,"slod"])
	m_mlod = max(ccQTL$qtls_alt[,"mlod"])

	# lod plots
	pdf("/g/steinmetz/brooks/genphen/growth/plots/ccQTL_lods.pdf")
			par(mfrow=c(2,1))
			plot(ccQTL$qtls_alt, ylim=c(0,m_slod), main=paste("ccQTL SLOD"),bandcol="gray90")
			abline(h=summary(ccQTL$permout)["5%","slod"], col="red", lty=2)
			abline(h=summary(ccQTL$permout)["10%","slod"], col="blue", lty=3)
			legend("topright", y.leg[i], c("5% FDR","10% FDR"), lty = c(2, 3), col = c("red","blue"))
      plot(ccQTL$qtls_alt, ylim=c(0,m_mlod), main=paste("ccQTL MLOD"),bandcol="gray90")
			abline(h=summary(ccQTL$permout)["5%","mlod"], col="red", lty=2)
			abline(h=summary(ccQTL$permout)["10%","mlod"], col="blue", lty=3)
			legend("topright", y.leg[i], c("5% FDR","10% FDR"), lty = c(2, 3), col = c("red","blue"))
	dev.off()

} else {
	#load(f)
}

#-------------------------------------------------------------------#
# growthrate with funQTL - MLOD/SLOD
# growth rate by regression using lm
#
#-------------------------------------------------------------------#
# rate from linear model
gr = filter(endometabolite, metabolite == "AKG", time_format == "absolute" ) %>%
     select(., strain, replicate, time, cellconc) %>%
     group_by(., strain) %>%
     do({
       x = as.numeric(.$time)
       y = as.numeric(log2(.$cellconc))
       o = lm(y~x)
       return(data.frame(rate = as.numeric(o$coefficients[2])))
       })
gr = data.frame(rate = gr$rate, row.names = sapply(as.character(gr$strain),strainRename))

f = "/g/steinmetz/brooks/genphen/growth/qtls/grQTL.rda"

if (FALSE) {
#if (!file.exists(f)) {
	# group phenotype timepoints

	cross_f = "/g/steinmetz/brooks/genphen/growth/cross_template.rda"
	#cross_f_abs = "/g/steinmetz/brooks/genphen/metabolome/cross_template_abs.rda"
	if (!file.exists(cross_f)) {
		cross =	runQTL(
					genotype = geno,
					phenotype = gr,
					marker_info = mrk,
					return_cross = TRUE,
					estimate.map = FALSE
					)
		save(cross, file = cross_f)
	} else {
		load(cross_f)
	}

	cross = calc.genoprob(cross, step = 0)
	# last phenotype column is the "id" tag
	pcols = seq(1, dim(cross$pheno)[2] - 1)
 	qtls = scanone(cross, method = "hk", pheno.col = pcols)
	eff = geteffects(cross, pheno.col = pcols)
	qtls_alt = scanoneF(cross, pheno.cols = pcols, method = "hk")
	# calc permutation threshold
	permout = scanoneF(cross, pheno.cols = pcols,
	                    method = "hk", n.perm = 1000, n.cluster = 20, verbose = F)
	# identify chrs with slod/lod above permute val
	chrs = unique(c(
		qtls_alt[qtls_alt[,"slod"]>=summary(permout)["5%","slod"],"chr"],
		qtls_alt[qtls_alt[,"mlod"]>=summary(permout)["5%","mlod"],"chr"]
		)
	)
	qtl_intervals = list()
	if (length(chrs)>0) {
		for (i in chrs) {
			qtl_intervals[[as.character(i)]] = mrk[rownames(bayesint(qtls_alt, chr = str_pad(i, 2, pad = "0"), prob=0.95, lodcolumn=1))]
		}
	}
	grQTL = list(qtls = qtls, eff = eff,
				qtls_alt = qtls_alt, qtl_intervals = qtl_intervals,
				permout = permout, pcols = pcols)
	save(grQTL, file = f)
	# plot stuff

	# effect plots
	pdf("/g/steinmetz/brooks/genphen/growth/plots/grQTL_effects.pdf")
		plotlod(grQTL$qtls, grQTL$eff, grQTL$pcols, gap=25, ylab="Time")
		mtext("grQTL", side = 3)
	dev.off()

	m_slod = max(grQTL$qtls_alt[,"slod"])
	m_mlod = max(grQTL$qtls_alt[,"mlod"])

	# lod plots
	pdf("/g/steinmetz/brooks/genphen/growth/plots/grQTL_lods.pdf")
			par(mfrow=c(2,1))
			plot(grQTL$qtls_alt, ylim=c(0,m_slod), main=paste("grQTL SLOD"),bandcol="gray90")
			abline(h=summary(grQTL$permout)["5%","slod"], col="red", lty=2)
			abline(h=summary(grQTL$permout)["10%","slod"], col="blue", lty=3)
			legend("topright", y.leg[i], c("5% FDR","10% FDR"), lty = c(2, 3), col = c("red","blue"))
      plot(grQTL$qtls_alt, ylim=c(0,m_mlod), main=paste("grQTL MLOD"),bandcol="gray90")
			abline(h=summary(grQTL$permout)["5%","mlod"], col="red", lty=2)
			abline(h=summary(grQTL$permout)["10%","mlod"], col="blue", lty=3)
			legend("topright", y.leg[i], c("5% FDR","10% FDR"), lty = c(2, 3), col = c("red","blue"))
	dev.off()

} else {
	#load(f)
}
