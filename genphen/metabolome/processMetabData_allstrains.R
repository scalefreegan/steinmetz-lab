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
library(plyr)
library(dplyr)
library(reshape2)
library(LSD)
library(qtl)
library(pheatmap)
library(funqtl)
library(parallel)
options(mc.cores = 24)
library(snow)
library(igraph)
#library(clustQTL)
#devtools::install_github("scalefreegan/steinmetz-lab/clustQTL")

# source rQTL utilities
devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/QTL/rQTL.R")

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
cleanData = function(df, row.names = NULL, nap = 0.9) {
	# remove columns or rows with all NA values
	tokr = apply(df,1,function(i)sum(is.na(i))<=nap*length(i))
	# no colnames
	tokc = which(!is.na(colnames(df)))
	df = df[tokr,tokc]
	if (is.numeric(row.names)) {
		rownames(df) = df[,row.names]
		df = df[,seq(1,dim(df)[2])[-row.names]]
	}
	return(df)
}

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

processData = function( f1, f2, f3, startRow = 2,... ) {
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
	thisf3_cellconc = cleanData(read.xlsx(f3,sheetName="Cell concentration",
		startRow=startRow,header=T),row.names=1)
	colnames(thisf3_cellconc) = strainRename(colnames(thisf3_cellconc))
	thisf3_biovol = cleanData(read.xlsx(f3,sheetName="Bio volume",
		startRow=startRow,header=T),row.names=1)
	colnames(thisf3_biovol) = strainRename(colnames(thisf3_biovol))
	thisf3_singlecell = cleanData(read.xlsx(f3,sheetName="Single cell volume",
		startRow=startRow,header=T),row.names=1)
	colnames(thisf3_singlecell ) = strainRename(colnames(thisf3_singlecell ))
	pb <- txtProgressBar(min = 0, max = length(sheets), style = 3)
	o = do.call(rbind,lapply((1:length(sheets)),function(i){
		#print(i)
		setTxtProgressBar(pb, i)
		i = sheets[i]
		thisf1 = cleanData(read.xlsx(f1,sheetName=i,startRow=startRow,header=T))
		thisf2 = cleanData(read.xlsx(f2,sheetName=i,startRow=startRow,header=T))
		common = intersect(colnames(thisf1),colnames(thisf2))
		# keep strains only
		common = common[grep("^X",common)]
		to_r = do.call(rbind,lapply(common,function(j){
			#print(j)
			do.call(rbind,lapply(c("relative","absolute"),function(k){
				#print(k)
				j = strainRename(j)
				strain = strsplit(j,"_")[[1]][1]
				replicate = as.numeric(strsplit(j,"_")[[1]][2])
				metabolite = i
				if (k=="relative") {
					touse = thisf1
					# add additional phes - only annotated for abs measure
					g = FALSE
					cellconc = rep(NA,length(touse[,1]))
					biovol = rep(NA,length(touse[,1]))
					singlecellvol = rep(NA,length(touse[,1]))
				} else {
					touse = thisf2
					# add additional phes - only annotated for abs measure
					g = FALSE
					cellconc = thisf3_cellconc[,j]
					biovol = thisf3_biovol[,j]
					singlecellvol = thisf3_singlecell[,j]
				}
				time = touse[,1]
				time_format = k
				# there are several negative values / i just take the absolute value
				# check for factor
				value = touse[,names(j)]
				if (class(value)=="factor") {
					value = as.numeric(levels(value)[value])
				} else {
					value = as.numeric(value)
				}
				# negative values should be zero
				value[value<0] = 0
				if (is.na(value[1])) {
				#   time = time[2:length(time)]
				#   value = value[2:length(value)]
				# 	cellconc = cellconc[2:length(cellconc)]
				# 	biovol = biovol[2:length(biovol)]
				# 	singlecellvol = singlecellvol[2:length(singlecellvol)]
						normInd = 2
				} else {
					normInd = 1
				}
				# value.log2 = value
				# update 12/10
				value.log2 = log2(value + 1)
				relative.log2 = sapply(value.log2,function(i){i/value.log2[normInd]})
				derivative.log2 = c(NA,diff(value.log2)/diff(time))
				data.frame(strain,replicate,metabolite,time,time_format,value,
				           value.log2,relative.log2,derivative.log2,
									 cellconc,biovol,singlecellvol)
			}))
		}))
		# remove NA
		# to_r = to_r[!is.na(to_r$value),]
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
f3 = paste(data_dir, "Dynamic_Metabolome_growth_and_morphology_113_strains.xlsx", sep="")

endo_f = "/g/steinmetz/project/GenPhen/data/endometabolome/data/endometabolite_full_12102015.rda"
#endo_f = "~/Desktop/tmpdata/full_endometabolome/endometabolite_full_23082015.rda"

if (file.exists(endo_f)) {
  load(endo_f)
} else{
  endometabolite = processData(f1, f2, f3, startRow = 2)
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

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/replicates_perstrain_abstime.pdf")
	  par(mfrow = c(2,2))
		d = filter(endometabolite, time_format == "absolute") %>%
		select(replicate,time,strain,value.log2) %>%
		group_by(strain) %>%
		do({
			x = filter(.,replicate==1)
	    y = filter(.,replicate==2)
	    t = sort(intersect(x$time,y$time))
	    x.log2 = x$value.log2[x$time%in%t]
	    y.log2 = y$value.log2[y$time%in%t]
			if (length(x.log2)==length(y.log2)) {
				data.frame(x.log2 = x.log2, y.log2 = y.log2)
			} else {
				data.frame()
			}
			})
		ggplot(d, aes(x.log2, y.log2)) +
			geom_point() +
			facet_wrap(~ strain)


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

	# Plot cell size characteristics ---------------------------------------------------
	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/size_char_bystrain.pdf", width=11.5,height=8)
	d = filter(endometabolite, time_format == "absolute", metabolite == "AKG")
	par(mfrow = c(3,1))
	hist(d$cellconc,100,main="Cell Concentration: all time points, all strains", xlab = "[cells/ml] ")
	hist(d$biovol,100,main="Biovolume: all time points, all strains", xlab = "[ul/ml] ")
	hist(d$singlecellvol,100,main="Single cell volume: all time points, all strains", xlab = "[fl] ")
	dev.off()

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/size_char_rep.pdf", width=11.5,height=8)
	par(mfrow = c(3,1))
	hist(d$cellconc,100,main="Cell Concentration: all time points, all strains", xlab = "[cells/ml] ")
	hist(d$biovol,100,main="Biovolume: all time points, all strains", xlab = "[ul/ml] ")
	hist(d$singlecellvol,100,main="Single cell volume: all time points, all strains", xlab = "[fl] ")
	dev.off()

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/size_char_dynamics.pdf", width=11.5,height=8)

	ggplot(
		data = d,
		aes(x = time, y = cellconc)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ strain, scales="fixed")
	scale_x_continuous(name="Time (absolute)") +
	ylab(expression(paste("cells/ml"))) +
	ggtitle("Physical Characteristics: Cell Conc.")

	ggplot(
		data = d,
		aes(x = time, y = biovol)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ strain, scales="fixed")
	scale_x_continuous(name="Time (absolute)") +
	ylab(expression(paste("ul/ml"))) +
	ggtitle("Physical Characteristics: Biovolume.")

	ggplot(
		data = d,
		aes(x = time, y = singlecellvol)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ strain, scales="fixed")
	scale_x_continuous(name="Time (absolute)") +
	ylab(name=expression(paste("fl"))) +
	ggtitle("Physical Characteristics: Single cell volume")

	dev.off()

	pdf("/g/steinmetz/project/GenPhen/data/endometabolome/plots/size_char_dynamics_allstrains_cellconc.pdf", width=11.5,height=8)
	ggplot(
		data = d,
		aes(x = time, y = cellconc)) +
	geom_point() +
	geom_smooth() +
	scale_x_continuous(name="Time (relative)") +
	scale_y_continuous(name=expression(paste("[cells/ml]"))) +
	ylim(c(min(d$cellconc),max(d$cellconc))) +
 	facet_wrap( ~ strain, scales="free_y") +
	ggtitle("Physical Characteristics: Cell Conc.")
	dev.off()

}

#-------------------------------------------------------------------#
# exploring rep behavior, etc
#
#-------------------------------------------------------------------#
if (F) {
	rep_all = filter(endometabolite,time_format=="absolute") %>% group_by(strain) %>% do(
		{
			o1 = try(cor(.$value.log2[.$replicate==1],.$value.log2[.$replicate==2],use="pair"),silent=T)
			if (class(o1)=="try-error") {
				o1 = NA
			}
			o2 = try(sd(.$value.log2[.$replicate==1],.$value.log2[.$replicate==2],na.rm=T),silent=T)
			if (class(o2)=="try-error") {
				o2 = NA
			}
			return(data.frame(mCor = o1, mSD = o2))
		}
		)


}


#-------------------------------------------------------------------#
# detect QTLs the traditional way, rQTL
#
#-------------------------------------------------------------------#


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
  filter(.,time_format=="relative")
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

f = "/g/steinmetz/brooks/genphen/metabolome/qtls/mQTLs_combrep.rda"
#if (!file.exists(f)) {
if (FALSE) {
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
	#load(f)
	#mQTLs_combrep = qtl_list
	#rm("qtl_list")
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
if (FALSE) {
#if (!file.exists(f)) {
	# group phenotype timepoints

	rep2metabolites = sapply(rownames(pheno), function(i) {
		o = strsplit(i,"_")[[1]][1]
		return(o)
	})

	cross_f = "/g/steinmetz/brooks/genphen/metabolome/cross_template.rda"
	#cross_f_abs = "/g/steinmetz/brooks/genphen/metabolome/cross_template_abs.rda"
	if (!file.exists(cross_f)) {
		cross =	runQTL(
					genotype = geno,
					phenotype = t(pheno),
					marker_info = mrk,
					return_cross = TRUE,
					estimate.map = FALSE
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
		try({
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
					qtl_intervals[[as.character(i)]] = mrk[rownames(bayesint(qtls_alt, chr = i, prob=0.95, lodcolumn=1))]
				}
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
	})
	names(mQTLs_funqtl_2014)  = unique(rep2metabolites)
	save(mQTLs_funqtl_2014, file = f)
	close(pb)
	# plot stuff

	# effect plots
	pdf("/g/steinmetz/brooks/genphen/metabolome/plots/funqtl_2014/effects.pdf")
		for (i in 1:length(mQTLs_funqtl_2014)) {
			try({
				plotlod(mQTLs_funqtl_2014[[i]]$qtls, mQTLs_funqtl_2014[[i]]$eff, mQTLs_funqtl_2014[[i]]$pcols, gap=25, ylab="Time")
				mtext(names(mQTLs_funqtl_2014)[i], side = 3)
			})
		}
	dev.off()

	# effect plots: jpeg
	for (i in 1:length(mQTLs_funqtl_2014)) {
		try({
			jpeg(paste("/g/steinmetz/brooks/genphen/metabolome/plots/funqtl_2014/effects",names(mQTLs_funqtl_2014)[i],".jpeg",sep=""))
				plotlod(mQTLs_funqtl_2014[[i]]$qtls, mQTLs_funqtl_2014[[i]]$eff, mQTLs_funqtl_2014[[i]]$pcols, gap=25, ylab="Time")
				mtext(names(mQTLs_funqtl_2014)[i], side = 3)
			dev.off()
		})
	}

	m_slod = max(as.numeric(unlist(lapply(seq(1:length(mQTLs_funqtl_2014)),function(i){try(max(mQTLs_funqtl_2014[[i]]$qtls_alt[,'slod']))}))), na.rm = T)
	m_mlod = max(as.numeric(unlist(lapply(seq(1:length(mQTLs_funqtl_2014)),function(i){try(max(mQTLs_funqtl_2014[[i]]$qtls_alt[,'mlod']))}))), na.rm = T)

	# lod plots
	pdf("/g/steinmetz/brooks/genphen/metabolome/plots/funqtl_2014/lods.jpeg")
		for (i in 1:length(mQTLs_funqtl_2014)) {
			try({
				par(mfrow=c(2,1))
				plot(mQTLs_funqtl_2014[[i]]$qtls_alt, ylim=c(0,m_slod), main=paste(names(mQTLs_funqtl_2014)[i],"SLOD"),bandcol="gray90")
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["5%","slod"], col="red", lty=2)
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["10%","slod"], col="blue", lty=3)
				legend("topright", y.leg[i], c("5% FDR","10% FDR"), lty = c(2, 3), col = c("red","blue"))
				plot(mQTLs_funqtl_2014[[i]]$qtls_alt, lodcolumn=2, ylim=c(0,m_mlod),
						main=paste(names(mQTLs_funqtl_2014)[i],"MLOD"), bandcol="gray90")
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["5%","mlod"], col="red", lty=2)
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["10%","mlod"], col="blue", lty=3)
			})
		}
	dev.off()

	# lod plots: jpeg

	for (i in 1:length(mQTLs_funqtl_2014)) {
		jpeg(paste("/g/steinmetz/brooks/genphen/metabolome/plots/funqtl_2014/lods",names(mQTLs_funqtl_2014)[i],".jpeg",sep=""))
			try({
				par(mfrow=c(2,1))
				plot(mQTLs_funqtl_2014[[i]]$qtls_alt, ylim=c(0,m_slod), main=paste(names(mQTLs_funqtl_2014)[i],"SLOD"),bandcol="gray90")
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["5%","slod"], col="red", lty=2)
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["10%","slod"], col="blue", lty=3)
				legend("topright", y.leg[i], c("5% FDR","10% FDR"), lty = c(2, 3), col = c("red","blue"))
				plot(mQTLs_funqtl_2014[[i]]$qtls_alt, lodcolumn=2, ylim=c(0,m_mlod),
						main=paste(names(mQTLs_funqtl_2014)[i],"MLOD"), bandcol="gray90")
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["5%","mlod"], col="red", lty=2)
				abline(h=summary(mQTLs_funqtl_2014[[i]]$permout)["10%","mlod"], col="blue", lty=3)
			})
		dev.off()
	}
} else {
	load(f)
}

#-------------------------------------------------------------------#
# Stitch
#
#-------------------------------------------------------------------#
fs = "/g/steinmetz/brooks/genphen/resources/stitch.rda"
# Load ---------------------------------------------------
if (!file.exists(fs)) {
	p2c = as.data.frame(read.table("/g/steinmetz/brooks/genphen/resources/4932.protein_chemical.links.transfer.v4.0.tsv",
		sep = "\t", header=T))
	# remove 4932 header on proteins
	p2c$protein = sub("4932\\.","",p2c$protein)

	chemNames = as.data.frame(read.delim("/g/steinmetz/brooks/genphen/resources/chemical.aliases.v4.0.tsv",
		sep = "\t", header=T))

	# only keep chemNames in p2c
	chemNames = filter(chemNames,chemical %in% levels(p2c$chemical)[p2c$chemical])
	save(p2c, chemNames, file = fs)
} else {
	load(fs)
}

# Try to map metabolites to chemID ---------------------------------------------------

m = unique(endometabolite$metabolite)
# annotated by hand. ANB 9/11/2015
m = data.frame(name = m, alias = c(
	"CID000000051", #AKG
	"CID000643757", #CAN
	"CID000000311", #CIT
	"CID000005793", #GLUC
	"CID000444972", #FUM
	"CID000000525", #MAL
	"CID000001005", #PEP
	"CID000001060", #PYR
	"CID000001110", #SUC
	"CID000000602", #ALA
	"CID000006322", #ARG
	"CID000006267", #ASN
	"CID000000424", #ASP
	"CID000000594", #CYS
	"CID000005961", #GLN
	"CID000000611", #GLU
	"CID000000750", #GLY
	"CID000006274", #HIS
	"CID000000779", #HSE
	"CID000000791", #ILE
	"CID000006106", #LEU
	"CID000005962", #LYS
	"CID000006137", #MET
	"CID000000994", #PHE
	"CID000145742", #PRO
	"CID000005951", #SER
	"CID000006305", #TRP
	"CID000001153", #TYR
	"CID000006287"  #VAL
	))
rownames(m) = m$alias

# how many m are in p2c
# sum(m$alias%in%p2c$chemical)/length(m$alias)
# 100%

# only chemicals in genephen, only combined score, normalize by each chemical such that the sum of scores per chemical equals 1
genphen_stitch = filter(p2c, chemical%in%m$alias) %>% select(., chemical, protein, combined_score)
genphen_stitch$alias = m[levels(genphen_stitch$chemical)[genphen_stitch$chemical],"name"]
genphen_stitch = genphen_stitch %>% group_by(., chemical) %>% do({
	s = .$combined_score
	score = .$combined_score/sum(.$combined_score)
	return(data.frame(chemical = .$chemical, alias = .$alias, protein = .$protein, score = score))
	})
# > dim(genphen_stitch)
# [1] 8308   3
save(genphen_stitch, file="/g/steinmetz/brooks/genphen/resources/genphen_stitch.rda")

# make it a network
v_meta = data.frame(vertex = c(levels(genphen_stitch$alias), genphen_stitch$protein), size = 3)
v_meta$vtype =  levels(v_meta$vertex)[v_meta$vertex]%in%levels(genphen_stitch$alias)
v_meta$shape[v_meta$vtype==TRUE] = "square"
v_meta$shape[v_meta$vtype==FALSE] = "circle"
v_meta$color[v_meta$vtype==TRUE] = RColorBrewer::brewer.pal(3,"Pastel2")[2]
v_meta$color[v_meta$vtype==FALSE] = RColorBrewer::brewer.pal(3,"Pastel2")[1]
v_meta$size[v_meta$vtype==TRUE] = 6
v_meta$size[v_meta$vtype==FALSE] = 3
v_meta$label[v_meta$vtype==TRUE] = levels(v_meta$vertex[v_meta$vtype==TRUE])[v_meta$vertex[v_meta$vtype==TRUE]]
v_meta$label[v_meta$vtype==FALSE] = NA
v_meta = v_meta[!duplicated(v_meta),]
v_meta = v_meta[rev(seq(1:dim(v_meta)[1])),]
g = graph_from_data_frame(genphen_stitch%>%select(.,alias,protein,score,chemical), directed = FALSE, vertices = v_meta)

# Select sig QTLs ---------------------------------------------------
#
# Taken directly from mQTL_explorer
devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/mQTL_explorer/global.R")

input = list()
input$m = "AKG"
input$co = 5
input$bco = 95
# select chromosomes
chrs = unique(data[[input$m]]$qtl[data[[input$m]]$qtl[,type]>=summary(data[[input$m]]$permout[,type],input$co/100)[1],"chr"])
chrs = levels(chrs)[chrs]

lodcolumn = if(type=="mlod"){ 2 } else { 1 }
qtl_intervals = list()
if (length(chrs)>0) {
	for (i in chrs) {
		qtl_intervals[[i]] = try(mrk[rownames(bayesint(data[[input$m]]$qtl, chr = str_pad(i, 2, pad = "0"), prob=input$bci/100, lodcolumn=lodcolumn))],silent = T)
		if (class(qtl_intervals[[i]])=="try-error") {
			qtl_intervals[[i]] = NULL
		} else {
			nn = sapply(as.character(seqnames(qtl_intervals[[i]])),function(i){
				paste(substr(i,1,3),as.roman(substr(i,4,5)),sep="")
			})
			qtl_intervals[[i]] = renameSeqlevels(qtl_intervals[[i]],nn)
			qtl_intervals[[i]] = keepSeqlevels(qtl_intervals[[i]],unique(nn))
			qtl_intervals[[i]] = range(qtl_intervals[[i]])
			qtl_intervals[[i]] = as.data.frame(cdsByOverlaps(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,qtl_intervals[[i]], type = "any", columns = "gene_id"))
		}
	}
}
qtl_df = do.call(rbind,qtl_intervals)
if (length(qtl_df) != 0) {
	qtl_df$gene_id = unlist(qtl_df$gene_id)
	gname_t = unlist(gname[unlist(qtl_df$gene_id)])
	gname_t = data.frame(gene_id = names(gname_t), name = gname_t)
	dname_t = unlist(dname[unlist(qtl_df$gene_id)])
	dname_t = data.frame(gene_id = names(dname_t), alias = dname_t)
	dname_t_long = unlist(dname_long[unlist(qtl_df$gene_id)])
	dname_t_long = data.frame(gene_id = names(dname_t_long), desc = dname_t_long)
	qtl_df = merge(qtl_df,gname_t,by="gene_id",sort=F,all.x=T)
	qtl_df = merge(qtl_df,dname_t,by="gene_id",sort=F,all.x=T)
	qtl_df = merge(qtl_df,dname_t_long,by="gene_id",sort=F,all.x=T)
	qtl_df = qtl_df[,c("gene_id","name","seqnames","start","end","strand","alias","desc")]
	colnames(qtl_df) = c("Sys.Name","Name","Chr","Start","End","Strand","Alias","Desc")
	#rownames(qtl_df) = qtl_df[,"Sys.Name"]
	qtl_df = qtl_df[!duplicated(qtl_df),]

}

# Plot  ---------------------------------------------------
if (.plot) {
	# assses features of p2c
	pdf("/g/steinmetz/brooks/genphen/resources/plots/stitch_score_dist.pdf")
		m <- ggplot(p2c, aes(x=combined_score)) + geom_histogram()
		m
	dev.off()
	# per gene
	hc = p2c %>% group_by(.,chemical) %>% summarise(.,mean = mean(combined_score))
	colnames(hc)[1] = "name"
	hc$source = "chemical"
	gc = p2c %>% group_by(.,protein) %>% summarise(.,mean = mean(combined_score))
	colnames(gc)[1] = "name"
	gc$source = "protein"
	cc = rbind(hc,gc)

	pdf("/g/steinmetz/brooks/genphen/resources/plots/stitch_score_dist_bysub.pdf")
		m <- ggplot(cc, aes(x = mean, fill = source)) + geom_density(alpha = 0.2)
		m
	dev.off()

	# normalized weights, genphen network
	pdf("/g/steinmetz/brooks/genphen/resources/plots/stitch_genphen_normweights.pdf")
		m <- ggplot(genphen_stitch, aes(x=score)) + geom_histogram()
		m
	dev.off()
	pdf("/g/steinmetz/brooks/genphen/resources/plots/stitch_genphen_normweights_perchem.pdf")
		m <- ggplot(genphen_stitch, aes(x=score)) + geom_histogram() + facet_wrap(~ alias)
		m
	dev.off()

	# genphen stitch network
	pdf("/g/steinmetz/brooks/genphen/resources/plots/stitch_genphen_network.pdf")
		plot.igraph(g, vertex.size = V(g)$size, vertex.label = NA, vertex.shape = V(g)$shape, vertex.color = V(g)$color)
	dev.off()
}
