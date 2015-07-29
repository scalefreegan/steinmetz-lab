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
source("~/Documents/git/steinmetz-lab/genphen/metabnorm/normfunctions.R")

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

f1 = "~/Desktop/tmpdata/endometabolome/Endometabolome_1B_25A_sorted by cultivation phase.xlsx"
f2 = "~/Desktop/tmpdata/endometabolome/Endometabolome_25A_sorted by cultivation time.xlsx"

endo_f = "~/Desktop/tmpdata/endometabolome/endometabolite_14072015.rda"
if (file.exists(endo_f)) {
  load(endo_f)
} else{
  endometabolite = processData(f1, f2, startRow = 2)
  save(endometabolite,file="~/Desktop/tmpdata/endometabolome/endometabolite_14072015.rda")
}

calc_deriv = function(df) {
  # calc first order finite difference
  # order df by time
  df = arrange( df,time )
  endo_rate = diff( df$endo_quant_log_normalized )
  exo_rate = diff( df$exo_quant_log_normalized )
  time = sapply( seq( 2, length( df$time ) ), function(i) {
    # change to ceil to keep time labeling consistent
    ceiling( mean( as.numeric( levels( df$time )[ df$time[ c( i-1,i ) ] ] ) ) )
  } )
  out = data.frame( metabolite = df$metabolite[1], strain = df$strain[1], parent = df$parent[1], time = time, endo_rate = endo_rate, exo_rate = exo_rate )
  return(out)
}

rel1 = function(df) {
  # order df by time
  df = arrange( df,time )
  df$endo_quant_rel = df$endo_quant_log_normalized/filter(df,time==min(as.numeric(levels(df$time))))$endo_quant_log_normalized
  df$exo_quant_rel = df$exo_quant_log_normalized/filter(df,time==min(as.numeric(levels(df$time))))$exo_quant_log_normalized
  return(df)
}

# normalize using mixed model
# convert endometabolome data to matrix, use log2 data

# # inititalize matrix
# metabolites = unique(data$metabolite)
# samples = sort(unique(paste(data$strain,data$time,sep="_")))
# endo_matrix = matrix(0,nrow=length(metabolites),ncol=length(samples),dimnames=list(metabolites,samples))
#
# # fill matrices
# for (i in 1:dim(data)[1]) {
# 	s = paste(data[i,]$strain,data[i,]$time,sep="_")
# 	m = levels(data[i,]$metabolite)[data[i,]$metabolite]
# 	v1 = data[i,]$endo_quant_log
# 	v2 = data[i,]$exo_quant_log
# 	endo_matrix[m,s] = v1
# 	exo_matrix[m,s] = v2
# }
# # assuming that each batch represents a bio rep
# batch = factor(unlist(lapply(strsplit(colnames(endo_matrix),"_"),"[",1)))
# cohort = factor(sapply(unlist(lapply(strsplit(colnames(endo_matrix),"_"),"[",1)),strtrim,1))
#
# # remove metabolites with NAs
# endo_matrix = endo_matrix[-which(apply(apply(endo_matrix,2,is.na),1,sum)>0),]
# exo_matrix = exo_matrix[-which(apply(apply(exo_matrix,2,is.na),1,sum)>0),]
#
# normalized_endo = normalizeMixed(endo_matrix,cohort=cohort,batch=batch)
# normalized_exo = normalizeMixed(exo_matrix,cohort=cohort,batch=batch)

# # add normalized data back
# data = data %>%
# 	rowwise() %>%
# 	do({
# 		metabolite = levels(.$metabolite)[.$metabolite]
# 		strain = paste(.$strain,.$time,sep="_")
# 		if ( metabolite%in%rownames(normalized_endo) & strain%in%colnames(normalized_endo) ) {
# 			data.frame(.,endo_quant_log_normalized = normalized_endo[metabolite,strain])
# 		} else {
# 			data.frame(.,endo_quant_log_normalized = NA)
# 		}
# 	}) %>%
# 	do({
# 		metabolite = levels(.$metabolite)[.$metabolite]
# 		strain = paste(.$strain,.$time,sep="_")
# 		if ( metabolite%in%rownames(normalized_exo) & strain%in%colnames(normalized_exo) ) {
# 			data.frame(.,exo_quant_log_normalized = normalized_exo[metabolite,strain])
# 		} else {
# 			data.frame(.,exo_quant_log_normalized = NA)
# 		}
# 	})

# Replicate behavior ---------------------------------------------------

repdata = endometabolite %>%
  group_by(metabolite,strain) %>%
  filter(.,time_format=="relative") %>%
  do({
    x = filter(.,replicate==1)
    y = filter(.,replicate==2)
    t = sort(intersect(x$time,y$time))
    x = x$value.log2[x$time%in%t]
    y = y$value.log2[y$time%in%t]
    if (length(x)==length(y)) {
      data.frame(x=x,y=y)
    } else {
      data.frame()
    }
  })
pdf("~/Desktop/tmpdata/endometabolome/replicates_allmetabolites.pdf")
	heatscatter(repdata$x,repdata$y,main="Reproducibility, All Metabolites",ylab="Conc Log2(uM)",xlab="Conc Log2(uM)")
dev.off()
pdf("~/Desktop/tmpdata/endometabolome/replicates_permetabolite.pdf")
	par(mfrow = c(2,2))
	for (i in levels(repdata$metabolite)) {
		d = repdata[which(repdata$metabolite==i),]
		heatscatter(d$x,d$y,main=i,ylab="Conc Log2(uM)",xlab="Conc Log2(uM)",
			ylim=c(min(c(repdata$x,repdata$y),na.rm=T),max(c(repdata$x,repdata$y),na.rm=T)),
			xlim=c(min(c(repdata$x,repdata$y),na.rm=T),max(c(repdata$x,repdata$y),na.rm=T)))
	}
dev.off()

# Compute addtional stats ---------------------------------------------------

# calc relative levels
data = data %>%
	group_by( metabolite, strain ) %>%
	arrange() %>%
	do(rel1( . ) )

# calc derivative
data_dt = data %>%
	group_by( metabolite, strain ) %>%
	arrange() %>%
	do(calc_deriv( . ) )
data = merge( data, data_dt , all = TRUE )

# filter for metabolites measured in both strains
data_complete = data %>%
	group_by( metabolite ) %>%
	do({
	if( sum(!levels(.$parent) %in% .$parent)>=1) {
		df = data.frame()
	} else {
		#print("in")
		df = .
	}
})

endometabolite_long_relative = endometabolite%>%filter(.,time_format=="relative")%>%melt()

# Save data and env ---------------------------------------------------

save( data, data_complete, data_complete_long, file="./dynamic_metabolome_08052014/data.RData" )

# Plot overall data trends ---------------------------------------------------

pdf("~/Desktop/tmpdata/endometabolome/metabolome.pdf", width=11.5,height=8)


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



# everything expect s1/y1 data, make ref to s1 data for each column to plot reproducibility
rel21_data = filter( data, !strain=="S1" & !strain=="Y1" ) %>%
	rowwise() %>%
	do({
		parent = .$parent
		if (parent == "YJM789") {
			strain = "Y1"
		} else if (parent == "S288c") {
			strain = "S1"
		} else {
			print("FAIL!")
			return(data.frame())
		}
		filtered_data = filter(data, metabolite == .$metabolite, time == .$time)
		# don't know why i have to do this twice...
		filtered_data = filtered_data[filtered_data$strain==strain,]
		df = data.frame(., endo_ref = filtered_data$endo_quant, exo_ref = filtered_data$exo_quant, endo_ref_norm = filtered_data$endo_quant_log_normalized, exo_ref_norm = filtered_data$exo_quant_log_normalized)
		return(df)
		})

var_data = data %>%
	group_by(metabolite,parent,time) %>%
	summarise(
		endo_mean=mean(endo_quant_log_normalized,na.rm=TRUE),
		endo_variance=var(endo_quant_log_normalized,na.rm=TRUE),
		exo_mean=mean(exo_quant_log_normalized,na.rm=TRUE),
		exo_variance=var(exo_quant_log_normalized,na.rm=TRUE)
	)

pdf("./dynamic_metabolome_08052014/data_qc.pdf", width=11.5,height=8)
# endo
ggplot( rel21_data, aes( x = endo_ref, y = endo_quant, group = parent, colour = parent ) ) +
	geom_point() +
	geom_abline() +
    scale_y_continuous(name="Replicates, Metabolite Levels (uM)") +
    scale_x_continuous(name="S1, Metabolite Levels (uM)") +
    ggtitle("Endo-Metabolite Reproducibility")

ggplot( rel21_data, aes( x = log2( endo_ref ), y = log2( endo_quant ), group = parent, colour = parent ) ) +
	geom_point() +
	geom_abline() +
    scale_y_continuous(name="Replicates, Metabolite Levels log2(uM)") +
    scale_x_continuous(name="S1, Metabolite Levels log2(uM)")  +
    ggtitle("Endo-Metabolite Reproducibility, Log Scale")

ggplot( rel21_data, aes( x = endo_ref_norm, y = endo_quant_log_normalized, group = parent, colour = parent ) ) +
	geom_point() +
	geom_abline() +
    scale_y_continuous(name="Replicates, Metabolite Levels log2(uM)") +
    scale_x_continuous(name="S1, Metabolite Levels log2(uM)")  +
    ggtitle("Endo-Metabolite Reproducibility, Log Scale, Normalized")
ggplot( var_data, aes( x = endo_mean, y = endo_variance, group = parent, colour = parent ) ) +
	geom_point() +
	geom_smooth(method=glm, se=F) +
 	scale_y_continuous(name="Variance") +
    scale_x_continuous(name="Mean Metabolite Levels log2(uM), normalized")  +
    ggtitle("Endo-Metabolite Mean Dependent Variance, Log Scale, Normalized")
# exo
ggplot( rel21_data, aes( x = exo_ref, y = exo_quant, group = parent, colour = parent ) ) +
	geom_point() +
	geom_abline() +
    scale_y_continuous(name="Replicates, Metabolite Levels (uM)") +
    scale_x_continuous(name="S1, Metabolite Levels (uM)") +
    ggtitle("Exo-Metabolite Reproducibility")

ggplot( rel21_data, aes( x = log2( exo_ref ), y = log2( exo_quant ), group = parent, colour = parent ) ) +
	geom_point() +
	geom_abline() +
    scale_y_continuous(name="Replicates, Metabolite Levels log2(uM)") +
    scale_x_continuous(name="S1, Metabolite Levels log2(uM)")  +
    ggtitle("Exo-Metabolite Reproducibility, Log Scale")

ggplot( rel21_data, aes( x = exo_ref_norm, y = exo_quant_log_normalized, group = parent, colour = parent ) ) +
	geom_point() +
	geom_abline() +
    scale_y_continuous(name="Replicates, Metabolite Levels log2(uM)") +
    scale_x_continuous(name="S1, Metabolite Levels log2(uM)")  +
    ggtitle("Exo-Metabolite Reproducibility, Log Scale, Normalized")

ggplot( var_data, aes( x = exo_mean, y = exo_variance, group = parent, colour = parent ) ) +
	geom_point() +
	geom_smooth(method=glm, se=F) +
 	scale_y_continuous(name="Variance") +
    scale_x_continuous(name="Mean Metabolite Levels log2(uM), normalized")  +
    ggtitle("Exo-Metabolite Mean Dependent Variance, Log Scale, Normalized")
dev.off()



