#! /usr/bin/env Rscript
# designed to be saved as .Rprofile in working dir
# devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/genphen/metabolome/processMetabData_allstrains.R")
#-------------------------------------------------------------------#
# Process Metabolome Data Nicole on 07.07.2015
# ! Transform into a computable format !
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
library(xlsx)
library(ggplot2)
#library(plotly)
library(dplyr)
library(reshape2)
library(LSD)
library(qtl)
library(pheatmap)
library(funqtl)
library(parallel)
library(snow)
#library(clustQTL)
#devtools::install_github("scalefreegan/steinmetz-lab/clustQTL")

# Import functions ---------------------------------------------------
# not required currently
# source("/g/steinmetz/brooks/git/steinmetz-lab/genphen/metabnorm/normfunctions.R")

# NOTE: this function is specific to data format supplied by Nicole, which looks like this:
# Intracellular concentration [µM]
# Cultivation time [h]	1B_1	1B_2	1C_1	1C_2	1D_1	1D_2
# 16
# 17	657.4262347	767.9153504	389.2626001	405.617732	684.0041913	817.4367778
# 18	955.6472403	1118.252453	342.7908314	363.6082801	660.4022457	671.3698753
# 19	1733.025721	1892.588403	428.7290925	536.6023722	652.8111875	680.7000702
# 20	3098.518239	3541.619307	1203.66303	1277.799992	602.1646668	481.5027229
# This is slight different than the original format. Each of the metabolites is contained in a separate worksheet

# Custom functions ---------------------------------------------------

processData = function( f1, f2, startRow = 2,... ) {
	# f1 is relative time - defined by "cultivation phase"
	# eg. Endometabolome_1B_25A_sorted by cultivation phase.xlsx
	# f2 contains measurements at absolute time
	# eg. Endometabolome_25A_sorted by cultivation time.xlsx
	# f3 contains additional measurements for the strains
	# including Cell_conc, biovolume, single cell volume, and traits
	# we will not use "traits" sheet - will calculate them myself if needed
	wb1 = loadWorkbook(f1)
	wb2 = loadWorkbook(f2)
	wb3 = loadWorkbook(f3)
	sheets = intersect(names(getSheets(wb1)),names(getSheets(wb2)))
	sheets.g = names(getSheets(wb3))[names(getSheets(wb3))!="Traits"]
	pb <- txtProgressBar(min = 0, max = length(sheets), style = 3)
	o = do.call(rbind,lapply((1:length(sheets)),function(i){
		#print(i)
		setTxtProgressBar(pb, i)
		i = sheets[i]
		thisf1 = read.xlsx(f1,sheetName=i,startRow=startRow,header=T)
		thisf2 = read.xlsx(f2,sheetName=i,startRow=startRow,header=T)
		thisf3 = read.xlsx(f2,sheetName=i,startRow=startRow,header=T)
		common = intersect(colnames(thisf1),colnames(thisf2))
		# keep strains only
		common = common[grep("^X",common)]
		to_r = do.call(rbind,lapply(common,function(j){
			do.call(rbind,lapply(c("relative","absolute"),function(k){
				strain = strsplit(j,"_")[[1]][1]
				if (nchar(strain)==3) {
					strain = gsub("^X","0",strain)
				} else {
					strain = gsub("^X","",strain)
				}
				replicate = as.numeric(strsplit(j,"_")[[1]][2])
				metabolite = i
				if (k=="relative") {
					touse = thisf1
					# add additional phes - only annotated for abs measure
					g = FALSE
				} else {
					touse = thisf2
					# add additional phes - only annotated for abs measure
					g = FALSE
				}
				time = touse[,1]
				time_format = k
				value = abs(as.numeric(touse[,j]))
				if (is.na(value[1])) {
				  time = time[2:length(time)]
				  value = value[2:length(value)]
				}
				value.log2 = value
				value.log2[value>0 & !is.na(value)] = log2(value.log2)
				relative.log2 = sapply(value.log2,function(i){i/value.log2[1]})
				derivative.log2 = c(NA,diff(value.log2)/diff(time))
				data.frame(strain,replicate,metabolite,time,time_format,value,
				           value.log2,relative.log2,derivative.log2)
			}))
		}))
		# remove NA
		to_r = to_r[!is.na(to_r$value),]
		return(to_r)
	}))
	close(pb)
	return(o)
}

# Load and process data ---------------------------------------------------

#data_dir = "~/Desktop/tmpdata/full_endometabolome/"
data_dir = "/g/steinmetz/project/GenPhen/data/endometabolome/data/"

f1 = paste(data_dir, "Endometabolome_1B_46B_sorted by cultivation phase.xlsx", sep="")
f2 = paste(data_dir, "Endometabolome_46B_sorted by cultivation time.xlsx", sep="")
f3 = paste(data_dir, "Dynamic Metabolome_growth and morphology_113 strains.xlsx", sep="")

endo_f = "/g/steinmetz/project/GenPhen/data/endometabolome/endometabolite_full_23082015.rda"
#endo_f = "~/Desktop/tmpdata/full_endometabolome/endometabolite_full_23082015.rda"

if (file.exists(endo_f)) {
  load(endo_f)
} else{
  endometabolite = processData(f1, f2, startRow = 2)
	# fix replicate problem: change 28B_1/28B_1 to 28B_1/28B_2
	endometabolite[which(endometabolite$replicate==1.1),"replicate"] = 2
	endometabolite$replicate = as.factor(endometabolite$replicate)
	# replace NaN with NA

  save(endometabolite, file = endo_f)
}

# Replicate behavior ---------------------------------------------------
if (.plot) {
	repdata = endometabolite %>%
	  group_by(metabolite,strain) %>%
	  filter(.,time_format=="relative") %>%
	  do({
	    x = filter(.,replicate==1)
	    y = filter(.,replicate==2)
	    t = sort(intersect(x$time,y$time))
	    x.log2 = x$value.log2[x$time%in%t]
	    y.log2 = y$value.log2[y$time%in%t]
			x.diffBYmean = x$derivative.log2[x$time%in%t]/mean(x.log2)
	    y.diffBYmean = y$derivative.log2[y$time%in%t]/mean(y.log2)
	    if (length(x.log2)==length(y.log2)) {
	      data.frame(x.log2 = x.log2, y.log2 = y.log2,
					x.diffBYmean = x.diffBYmean, y.diffBYmean = y.diffBYmean)
	    } else {
	      data.frame()
	    }
	  })

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/replicates_allmetabolites.pdf")
		heatscatter(repdata$x.log2,repdata$y.log2,main="Reproducibility,
		All Metabolites, Paired Time",
			ylab="[X]] Log2(uM), Rep 2",xlab="[X] Log2(uM), Rep 1")
		# remove outliers
		lim1 = max(quantile(repdata$x.diffBYmean,na.rm=T,0.99),
			quantile(repdata$y.diffBYmean,na.rm=T,0.99))
		lim2 = max(quantile(repdata$x.diffBYmean,na.rm=T,0.01),
				quantile(repdata$y.diffBYmean,na.rm=T,0.01))
		d = repdata %>% filter(., x.diffBYmean <= lim1 & x.diffBYmean >= lim2 &
			y.diffBYmean <= lim1 & y.diffBYmean >= lim2)
		heatscatter(d$x.diffBYmean,d$y.diffBYmean,
			main="Reproducibility, All Metabolites, Paired Time, Derivative By Mean",
			ylab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 2",
			xlab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 1"
		)
	dev.off()

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/replicates_permetabolite.pdf")
		par(mfrow = c(2,2))
		for (i in levels(repdata$metabolite)) {
			d = repdata[which(repdata$metabolite==i),]
			heatscatter(d$x.log2,d$y.log2,main=i,ylab="[X] Log2(uM), Rep 2",xlab="[X] Log2(uM) Rep 1",
				ylim=c(min(c(repdata$x.log2,repdata$y.log),na.rm=T),
					max(c(repdata$x.log2,repdata$y.log),na.rm=T)),
				xlim=c(min(c(repdata$x.log2,repdata$y.log),na.rm=T),
					max(c(repdata$x.log2,repdata$y.log),na.rm=T))
			)
			d = repdata[which(repdata$metabolite==i),]
			lims = max(quantile(repdata$x.diffBYmean,na.rm=T,0.99),
				quantile(repdata$y.diffBYmean,na.rm=T,0.99))
			heatscatter(d$x.diffBYmean,d$y.diffBYmean,
				main = i,
				ylab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 2",
				xlab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 1",
				xlim = c(-lims,lims),
				ylim = c(-lims,lims)
			)
		}
	dev.off()

	# Plot overall data trends ---------------------------------------------------

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/metabolome.pdf", width=11.5,height=8)
	ggplot(
		data = endometabolite %>%
			filter(! is.na(value.log2)) %>%
			filter(time_format=="relative") %>%
			group_by(metabolite) %>%
			do({filter(.,abs(.$value.log2-mean(.$value.log2))<5*sd(.$value.log2))}),
		aes(x = time, y = value.log2)) +
	  	geom_point() +
		geom_smooth() +
		facet_wrap( ~ metabolite, scales="free_y") +
		scale_x_discrete(name="Time (relative)") +
	    scale_y_continuous(name=expression(paste("Level log2(",mu,"M)"))) +
	    ggtitle("Endo- Metabolome Levels Across All Strains All Metabolites")
	dev.off()

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/metabolome_genotypefreq.pdf", width=11.5,height=8)
		hist(apply(geno,1,function(i)sum(i==1)/length(i)),main="Genotype composition, all markers",xlab="Genotype == 1, Freq",xlim=c(0,1))
	dev.off()

	repdata_abs = endometabolite %>%
	  group_by(metabolite,strain) %>%
	  filter(.,time_format=="absolute") %>%
	  do({
	    x = filter(.,replicate==1)
	    y = filter(.,replicate==2)
	    t = sort(intersect(x$time,y$time))
	    x.log2 = x$value.log2[x$time%in%t]
	    y.log2 = y$value.log2[y$time%in%t]
			x.diffBYmean = x$derivative.log2[x$time%in%t]/mean(x.log2)
	    y.diffBYmean = y$derivative.log2[y$time%in%t]/mean(y.log2)
	    if (length(x.log2)==length(y.log2)) {
	      data.frame(x.log2 = x.log2, y.log2 = y.log2,
					x.diffBYmean = x.diffBYmean, y.diffBYmean = y.diffBYmean)
	    } else {
	      data.frame()
	    }
	  })

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/replicates_allmetabolites_abstime.pdf")
	  heatscatter(repdata_abs$x.log2,repdata_abs$y.log2,main="Reproducibility,
	  All Metabolites, Paired Time",
	    ylab="[X]] Log2(uM), Rep 2",xlab="[X] Log2(uM), Rep 1")
	  # remove outliers
	  lim1 = max(quantile(repdata_abs$x.diffBYmean,na.rm=T,0.99),
	    quantile(repdata_abs$y.diffBYmean,na.rm=T,0.99))
	  lim2 = max(quantile(repdata_abs$x.diffBYmean,na.rm=T,0.01),
	      quantile(repdata_abs$y.diffBYmean,na.rm=T,0.01))
	  d = repdata_abs %>% filter(., x.diffBYmean <= lim1 & x.diffBYmean >= lim2 &
	    y.diffBYmean <= lim1 & y.diffBYmean >= lim2)
	  heatscatter(d$x.diffBYmean,d$y.diffBYmean,
	    main="Reproducibility, All Metabolites, Paired Time, Derivative By Mean",
	    ylab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 2",
	    xlab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 1"
	  )
	dev.off()

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/replicates_permetabolite_abstime.pdf")
	  par(mfrow = c(2,2))
	  for (i in levels(repdata_abs$metabolite)) {
	    d = repdata_abs[which(repdata_abs$metabolite==i),]
	    heatscatter(d$x.log2,d$y.log2,main=i,ylab="[X] Log2(uM), Rep 2",xlab="[X] Log2(uM) Rep 1",
	      ylim=c(min(c(repdata_abs$x.log2,repdata_abs$y.log),na.rm=T),
	        max(c(repdata_abs$x.log2,repdata_abs$y.log),na.rm=T)),
	      xlim=c(min(c(repdata_abs$x.log2,repdata_abs$y.log),na.rm=T),
	        max(c(repdata_abs$x.log2,repdata_abs$y.log),na.rm=T))
	    )
	    d = repdata_abs[which(repdata_abs$metabolite==i),]
	    lims = max(quantile(repdata_abs$x.diffBYmean,na.rm=T,0.99),
	      quantile(repdata_abs$y.diffBYmean,na.rm=T,0.99))
	    heatscatter(d$x.diffBYmean,d$y.diffBYmean,
	      main = i,
	      ylab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 2",
	      xlab = "(d[X]/dt)/mean([X]), Log2(uM), Rep 1",
	      xlim = c(-lims,lims),
	      ylim = c(-lims,lims)
	    )
	  }
	dev.off()

	# Plot overall data trends ---------------------------------------------------

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/metabolome_abstime.pdf", width=11.5,height=8)
	ggplot(
	  data = endometabolite %>%
	    filter(! is.na(value.log2)) %>%
	    filter(time_format=="relative") %>%
	    group_by(metabolite) %>%
	    do({filter(.,abs(.$value.log2-mean(.$value.log2))<5*sd(.$value.log2))}),
	  aes(x = time, y = value.log2)) +
	    geom_point() +
	  geom_smooth() +
	  facet_wrap( ~ metabolite, scales="free_y") +
	  scale_x_discrete(name="Time (relative)") +
	    scale_y_continuous(name=expression(paste("Level log2(",mu,"M)"))) +
	    ggtitle("Endo- Metabolome Levels Across All Strains All Metabolites")
	dev.off()

}


#-------------------------------------------------------------------#
# detect QTLs the traditional way, rQTL
#
#-------------------------------------------------------------------#

# source rQTL utilities
devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/QTL/rQTL.R")
#
# Don't think this is necessary for rqtl package
# in fact, I think it will cause imputaition of genotypes, which I would like to
# avoid
#
# genotype_f = "/g/steinmetz/project/GenPhen/data/endometabolome/genotypesANDmarkers_allmetabolites_23082015.rda"
# if (file.exists(genotype_f)) {
#   load(genotype_f)
# } else{
# 	load("/g/steinmetz/brooks/yeast/genomes/S288CxYJM789/genotypes_S288c_R64.rda")
# 	strains = levels(endometabolite$strain)
# 	gm_list = clustQTL::combineMarkers(geno, mrk, limit_strains = strains, limit_markers = NULL,
# 	                          na.rm = TRUE, rm_type = c("marker","strain")[1], collpase_markers = TRUE,
# 	                          marker_rename = FALSE, impute_markers = TRUE, clean_markers = TRUE)
# 	geno = gm_list$genotypes
# 	mrk = gm_list$markers
#   save(geno,mrk, file = genotype_f)
# }

# load genotype and markers files
genotype_f = "/g/steinmetz/brooks/yeast/genomes/S288CxYJM789/genotypes_S288c_R64.rda"
load(genotype_f)

# 7/10/2015 changed from time_format=="relative" to time_format=="absolute"
data = endometabolite %>%
  #filter(metabolite=="ARG") %>%
  filter(.,time_format=="absolute")
pnames = rep(NA,length(unique(paste(data$metabolite,data$replicate,data$time,sep="_"))))
names(pnames) = unique(paste(data$metabolite,data$replicate,data$time,sep="_"))

repnames = unique(paste(data$metabolite,data$time,sep="_"))
nreps = length(unique(data$replicate))

# with reps as separate measurements
# This way detects more QTLs but I'm 95% sure it is
# statistically problematic
# pheno = do.call(cbind,lapply(levels(data$strain),function(i){
# 		strain_data = filter(data,strain==i)
# 		# combine replicates
# 		vals = do.call(cbind,lapply(repnames,function(j){
# 			#print(j)
# 			to_r = data.frame(value.log2 = rep(NA,nreps))
# 			repj = strsplit(j,split="_")[[1]]
# 			m = repj[1]
# 			tj = repj[2]
# 			thisv = filter(strain_data,metabolite == m,time == tj)[,"value.log2"]
# 			if (length(thisv) > 0) {
# 					to_r$value.log2[seq(1,length(thisv))] = thisv
# 			}
# 			return(to_r)
# 			}))
# 		vals = t(vals)
# 		rownames(vals) = repnames
# 		colnames(vals) = rep(i,dim(vals)[2])
# 		return(vals)
# 	}))

# take mean of replicates
pheno = do.call(cbind,lapply(levels(data$strain),function(i){
		strain_data = filter(data,strain == i)
		# combine replicates
		vals = sapply(repnames,function(j){
			#print(j)
			mj = strsplit(j,"_")[[1]][1]
			tj = as.numeric(strsplit(j,"_")[[1]][2])
			mean(
				filter(data, strain == i, metabolite == mj, time == tj)[,"value.log2"],
				na.rm =T)
			})
		names(vals) = repnames
		return(vals)
	}))
colnames(pheno) = levels(data$strain)

# change NaN to NA
pheno[is.nan(pheno)] = NA

# order rownames by time
pheno = pheno[order(rownames(pheno)),]

# remove missing data
# pheno = pheno[which(apply(pheno,1,function(i)sum(is.na(i)))==0),]

# combine all data per metabolite
mnames = levels(data$metabolite)
pheno_comb = do.call(cbind,lapply(mnames,function(i){
	i_ind = which(sapply(rownames(pheno), function(j){strsplit(j,"_")[[1]][1]})%in%i)
	do.call(rbind,lapply(seq(1,dim(pheno[i_ind,])[1]),function(j) t(pheno[i_ind,][j,,drop=F])))
}))
colnames(pheno_comb) = sapply(colnames(pheno_comb),function(i){strsplit(i,"_")[[1]][1]})

f = "/g/steinmetz/brooks/genphen/metabolome/mQTLs_combrep.rda"
if (!file.exists(f)) {
	mQTLs_combrep =	runQTL(
		genotype = geno,
    phenotype = t(pheno),
    marker_info = mrk,
    permute = T, # compute significance of each QTL LOD by permutation
    pca = F, # maximize QTL detection by removing confounders using PCA
    permute_alpha = 0.05,
    save_file = f)
	write.table(do.call(rbind, mQTLs_combrep$sig_qtls),
		file = gsub(".rda","_sigtable.txt",f),
		sep = "\t",
		row.names = F,
		col.names = T
		)

	# plot lod profiles for each metabolite
	pdf(gsub(".rda","_lod_profiles.pdf",f), width=11.5,height=8)
	png(gsub(".rda","_lod_profiles.png",f), width=20,height=8)
		for (i in 1:length(mQTLs_combrep$qtls)) {
			x = mQTLs_combrep$qtls[[i]]
			mx = max(x[,3])
			if (mx < 4) {
				mx = 4
			}
			x[,3] = 10^-x[,3]
			colnames(x)[3] = "pval"
			plotManhattan(qtls = x, mrk  = mrk, main = names(mQTLs_combrep$qtls)[i], ylim = c(0,mx),suggestiveline = mQTLs_combrep$qtls_threshold[[i]])
		}
	dev.off()

	# plot correlations between genome wide lod scores for all metabolite profiles
	# correlation of LOD scores between timepoints
	cor_m = do.call(rbind, lapply(seq(1:length(mQTLs_combrep$qtls)),function(i){
			sapply(seq(1:length(mQTLs_combrep$qtls)),function(j){
					cor(mQTLs_combrep$qtls[[i]][,"lod"], mQTLs_combrep$qtls[[j]][,"lod"])
				})
		}))
	rownames(cor_m) = names(mQTLs_combrep$qtls)
	colnames(cor_m) = names(mQTLs_combrep$qtls)
	pdf(gsub(".rda","_lod_timepointcor.pdf",f), width=11.5,height=8)
		pheatmap(cor_m,breaks=seq(-1, 1, length.out = 100), fontsize = 6, main = "Correlation LOD score, Metabolites/Timepoints" )
	dev.off()

	# plot pheno dist for each sig qtl
	pdf(gsub(".rda","_sigqtl_phenos.pdf",f), width=11.5,height=8)
		sigqtls = do.call(rbind, mQTLs_combrep$sig_qtls)
		to_p = lapply(1:dim(sigqtls)[1], function(i) {
			m = sigqtls[i,"i"]
			n = sigqtls[i,"names"]
			g = geno[n,]
			x = pheno[m,]
			df = data.frame(geno = g, x = x)
			p <- ggplot(df, aes(factor(geno), x, fill = factor(geno))) +
				geom_boxplot(na.rm=T) +
				geom_jitter(position = position_jitter(width = .1), na.rm=T) +
				scale_x_discrete(name="Genotype") +
				scale_y_continuous(name=expression("[X] Log2(uM))")) +
	    	ggtitle(paste(m,n))
			p
		})
		to_p[1:length(to_p)]
	dev.off()
} else {
	load(f)
	mQTLs_combrep = qtl_list
	rm("qtl_list")
}

#-------------------------------------------------------------------#
# Archive!
#
#-------------------------------------------------------------------#
# # NOTE: this is probably statistically invalid way to combine across markers
# # do not use it for any serious business
# f = "/g/steinmetz/brooks/genphen/dynamic_metabolome_20082015/mQTLs_comball.rda"
# if (!file.exists(f)) {
# 	mQTLs_permetabolite =	runQTL(
# 		genotype = geno,
#     phenotype = pheno_comb,
#     marker_info = mrk,
#     permute = T, # compute significance of each QTL LOD by permutation
#     pca = F, # maximize QTL detection by removing confounders using PCA
#     permute_alpha = 0.05,
#     save_file = f)
# 	write.table(do.call(rbind, mQTLs_permetabolite$sig_qtls),
# 		file = gsub(".rda","_sigtable.txt",f),
# 		sep = "\t",
# 		row.names = F,
# 		col.names = T
# 		)
# } else {
# 	#load(f)
# }

#-------------------------------------------------------------------#
# With funQTL - MLOD/SLOD
# No smoothing / PCA
#
#-------------------------------------------------------------------#
f = "/g/steinmetz/brooks/genphen/metabolome/qtls/mQTLs_comball_funqtl_2014.rda"
if (F) {
#if (!file.exists(f)) {
	# group phenotype timepoints

	rep2metabolites = sapply(rownames(pheno), function(i) {
		o = strsplit(i,"_")[[1]][1]
		return(o)
	})

	cross_f = "/g/steinmetz/brooks/genphen/metabolome/cross_template.rda"
	if (!file.exists(cross_f)) {
		cross =	runQTL(
					genotype = geno,
					phenotype = t(pheno),
					marker_info = mrk,
					return_cross = TRUE,
					estimate.map=FALSE
					)
		save(cross, file = cross_f)
	} else {
		load(cross_f)
	}

	# first detect qtls individually
	# have to run funqtl for each metabolite
	# as shortcut, just replace cross object phenotypes
	# for each metabolite
	pb = txtProgressBar(min = 0, max = length(unique(rep2metabolites)), style = 3)
	mQTLs_funqtl_2014 = mclapply(1:length(unique(rep2metabolites)), function(i) {
		m = unique(rep2metabolites)[i]
		setTxtProgressBar(pb, i)
		these_phe = c(names(which(rep2metabolites == m)),"id")
		cross_tmp = cross
		cross_tmp$pheno = cross_tmp$pheno[,these_phe]
		colnames(cross_tmp$pheno) = c(sapply(colnames(cross_tmp$pheno), function(i){
			o = strsplit(i, split = "_")[[1]][2]
			if (is.na(o)) {
				o = "id"
			}
			return(o)
			}))
		cross_tmp = calc.genoprob(cross_tmp, step = 0)
		# last phenotype column is the "id" tag
		pcols = seq(1, length(these_phe) - 1)
	 	qtls = scanone(cross_tmp, method = "hk",pheno.col = pcols)
		eff = geteffects(cross_tmp, pheno.col = pcols)
		qtls_alt = scanoneF(cross_tmp, pheno.cols = pcols, method = "hk")
		# calc permutation threshold
		permout = scanoneF(cross_tmp, pheno.cols = pcols,
		                    method = "hk", n.perm = 1000, n.cluster = 20)
		# identify chrs with slod/lod above permute val
		chrs = unique(c(
			qtls_alt[qtls_alt[,"slod"]>=summary(permout)["5%","slod"],"chr"],
			qtls_alt[qtls_alt[,"mlod"]>=summary(permout)["5%","mlod"],"chr"]
			)
		)
		qtl_intervals = list()
		for (i in chrs) {
			qtl_intervals[[as.character(i)]] = mrk[rownames(bayesint(qtls_alt, chr = i, prob=0.95, lodcolumn=1))]
		}
		# for future!
		# qtlslod_slod = stepwiseqtlF(cross_tmp, pheno.cols = pcols,
		# 	max.qtl = 6, usec = "slod",
		# 	method = "hk",
		# 	penalties = c(2.02, 2.62, 1.74) )
		# qtls_refined_slod = refineqtlF(cross_tmp, usec = "slod"
		# 		pheno.cols = pcols, method = "hk", qtl = qtlslod, keeplodprofile = T)
		# qtlslod_mlod = stepwiseqtlF(cross_tmp, pheno.cols = pcols,
		# 	max.qtl = 6, usec = "mlod",
		# 	method = "hk",
		# 	penalties = c(2.02, 2.62, 1.74) )
		# qtls_refined_mlod = refineqtlF(cross_tmp, usec = "mlod"
		# 		pheno.cols = pcols, method = "hk", qtl = qtlslod, keeplodprofile = T)
		return(list(qtls = qtls, eff = eff,
			qtls_alt = qtls_alt, qtl_intervals = qtl_intervals,
			permout = permout, pcols = pcols))
	})
	names(mQTLs_funqtl_2014)  = unique(rep2metabolites)
	close(pb)
	# plot stuff

	# effect plots
	pdf("/g/steinmetz/brooks/genphen/metabolome/plots/effects_funqtl_2014.rda")
		for (i in 1:length(mQTLs_funqtl_2014)) {
			plotlod(mQTLs_funqtl_2014[[i]]$qtls, mQTLs_funqtl_2014[[i]]$eff, mQTLs_funqtl_2014[[i]]$pcols, gap=25, ylab="Time")
			mtext(names(mQTLs_funqtl_2014)[i], side = 3)
		}
	dev.off()

	pdf("/g/steinmetz/brooks/genphen/metabolome/plots/lods_funqtl_2014.rda")
		for (i in 1:length(mQTLs_funqtl_2014)) {
			par(mfrow=c(2,1))
			plot(mQTLs_funqtl_2014[[i]]$qtls_alt, ylim=c(0,3.5), main="SLOD",
					bandcol="gray90")
			abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["5%","slod"], col="red", lty=3)
			plot(mQTLs_funqtl_2014[[i]]$qtls_alt, lodcolumn=2, ylim=c(0,5),
					main="MLOD", bandcol="gray90")
			abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["5%","mlod"], col="red", lty=3)
		}
	dev.off()
} else {
	load(f)
}

#-------------------------------------------------------------------#
# Table output of sig qtls
#
#-------------------------------------------------------------------#




#-------------------------------------------------------------------#
# Plot summaries
#
#-------------------------------------------------------------------#
