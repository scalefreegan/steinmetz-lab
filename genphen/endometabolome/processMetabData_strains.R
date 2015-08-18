#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# Process Metabolome Data Nicole on 07.07.2015
# ! Transform into a computable format !
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
library(xlsx)
library(ggplot2)
library(plotly)
library(dplyr)
library(reshape2)
library(LSD)

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
	wb1 = loadWorkbook(f1)
	wb2 = loadWorkbook(f2)
	sheets = intersect(names(getSheets(wb1)),names(getSheets(wb2)))
	pb <- txtProgressBar(min = 0, max = length(sheets), style = 3)
	o = do.call(rbind,lapply(1:length(sheets),function(i){
		setTxtProgressBar(pb, i)
		i = sheets[i]
		thisf1 = read.xlsx(f1,sheetName=i,startRow=startRow,header=T)
		thisf2 = read.xlsx(f2,sheetName=i,startRow=startRow,header=T)
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
				} else {
					touse = thisf2
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

data_dir = "/g/steinmetz/project/GenPhen/data/endometabolome/data/"

f1 = paste(data_dir, "Endometabolome_1B_25A_sorted by cultivation phase.xlsx", sep="")
f2 = paste(data_dir, "Endometabolome_25A_sorted by cultivation time.xlsx", sep="")

endo_f = "/g/steinmetz/project/GenPhen/data/endometabolome/endometabolite_14072015.rda"
if (file.exists(endo_f)) {
  load(endo_f)
} else{
  endometabolite = processData(f1, f2, startRow = 2)
  save(endometabolite, file = endo_f)
}

# Replicate behavior ---------------------------------------------------

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



# Save data and env ---------------------------------------------------

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

# Plot data reproducibility ---------------------------------------------------

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

#-------------------------------------------------------------------#
# detect QTLs the traditional way, rQTL
#
#-------------------------------------------------------------------#

# source rQTL utilities
devtools::source_url("https://raw.githubusercontent.com/scalefreegan/steinmetz-lab/master/QTL/rQTL.R")

genotype_f = "/g/steinmetz/project/GenPhen/data/endometabolome/genotypesANDmarkers_14072015.rda"
if (file.exists(genotype_f)) {
  load(genotype_f)
} else{
	load("/g/steinmetz/brooks/yeast/genomes/S288CxYJM789/genotypes_S288c_R64.rda")
	strains = levels(endometabolite$strain)
	gm_list = clustQTL::combineMarkers(geno, mrk, limit_strains = strains, limit_markers = NULL,
	                          na.rm = TRUE, rm_type = c("marker","strain")[1], collpase_markers = TRUE,
	                          marker_rename = FALSE, impute_markers = TRUE, clean_markers = TRUE)
	geno = gm_list$genotypes
	mrk = gm_list$markers
  save(geno,mrk, file = genotype_f)
}

# load genotype and markers files
ARGdata = endometabolite %>%
  filter(metabolite=="ARG") %>%
  filter(.,time_format=="relative")
pnames = rep(NA,length(unique(paste(ARGdata$metabolite,ARGdata$replicate,ARGdata$time,sep="_"))))
names(pnames) = unique(paste(ARGdata$metabolite,ARGdata$replicate,ARGdata$time,sep="_"))

pheno = t(do.call(cbind,lapply(levels(ARGdata$strain),function(i){
		strain_data = filter(ARGdata,strain==i)
		vals = pnames
		vals[paste(strain_data$metabolite,strain_data$replicate,strain_data$time,sep="_")] = strain_data[,"value.log2"]
		return(vals)
	})))
rownames(pheno) = levels(ARGdata$strain)

# remove missing data
# pheno = pheno[which(apply(pheno,1,function(i)sum(is.na(i)))==0),]

mQTLs =	runQTL(
			genotype = geno,
	    phenotype = pheno,
	    marker_info = mrk,
	    permute = T, # compute significance of each QTL LOD by permutation
	    pca = F, # maximize QTL detection by removing confounders using PCA
	    permute_alpha = 0.05,
	    save_file = "")
