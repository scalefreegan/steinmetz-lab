#! /usr/bin/env Rscript

#-------------------------------------------------------------------#
# Process Metabolome Data Nicole
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
library( ggplot2 )
library( plotly )
library( dplyr )
library( reshape2 )

# Import functions ---------------------------------------------------
source( "/g/steinmetz/brooks/git/steinmetz-lab/genphen/othertools/metabnorm/normfunctions.R" )

# NOTE: this function is specific to data format supplied by Nicole, which looks like this:
# Intracellulare concentration [µM]
# 	ISOCIT											AKG
# 	S1	S2	S3	S4	S5	S6	S7	Y1	Y2	Y3	Y4	S1
# 1	3844.50	3584.21	3947.04	4127.05	4275.45	4423.63	3901.31	2153.74	1961.84	2237.32	2114.77	34.1
# 2	2575.07	3084.54	3259.4	3199.84	3026.13	3414.68	3010.6	2151.62	2310.91	2146.4	1957.09	49.3
# 3	2579.87	2098.17	2552.03	3352.06	2661.11	2612.89	2529.99	1581.68	1648.46	1807.51	1589.14	31.8
# 4	2153.68	2157.96	2081.88	2047.92	2089.03	2078.23	1968.13	1341.41	1277.06	1430.45	1237.51	79.68
# 5	1609.98	1692.87	1625.89	1511.57	1674.03	1644.89	1489.67	1198.9	1207.58	1238.84	1175.75	29.49

# Custom functions ---------------------------------------------------

processData = function( f, metaboliteRow = 1, strainRow = 2, dataRows = NULL, dataRowsLabels = NULL, ignore_col = 1, ... ) {
	data = read.table(f,sep="\t", ...)
	data = data[ , -ignore_col]
	# clean up
	data = unique( data )
	for ( i in dim( data )[1] ) {
		if ( sum(data[i,]=="") == dim(data)[2] ) {
			data = data[ -i, ]
		}
	}
	if ( is.null(dataRows) ) {
		dataRows = seq( dim(data)[1] )[ -c( metaboliteRow,strainRow ) ]
		if ( is.null(dataRowsLabels ) ) {
			dataRowsLabels = seq( 1, length( dataRows ) )
		}
		if ( length(dataRows) != length(dataRowsLabels) ) {
			cat("ERROR: Please provide equal number of data rows and labels")
			return(NULL)
		}
	}
	# fill in columns
	i = 1
	while( i<=dim( data )[2] ) {
		if(i == 1) {
			name = as.character( data[ metaboliteRow, i ] )
		} else {
			if( data[ metaboliteRow, i ] == "" ) {
				levels( data[ , i ] ) <- c( levels( data[ , i ] ), name )
				data[ metaboliteRow, i ] = name
			} else {
				name = as.character( data[ metaboliteRow, i ] )
			}
		}
		i = i + 1
	}
	# convert to long format data
	out = data.frame()
	for ( i in seq(1, dim( data )[2] ) ) {
		for ( j in seq(1, length( dataRows ) ) ) {
			if ( grepl("S",as.character( data[ strainRow,i ] ) ) ) {
				parent = "S288c"
			} else {
				parent = "YJM789"
			}
			out = rbind( out, data.frame(
				metabolite = as.character( data[ metaboliteRow,i ] ),
				strain =  as.character( data[ strainRow,i ] ),
				parent = parent,
				data_label = as.character( dataRowsLabels[ j ] ),
				data = as.numeric( levels( data[ dataRows[ j ],i ] ) )[ data[ dataRows[ j ],i ] ]
				) )
		}
	}
	# clean up
	out = out[ rowSums(is.na(out)) == 0, ]
	return(out)
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




# Load and process data ---------------------------------------------------

f1 = "./dynamic_metabolome_08052014/Endometabolome_ParentalStrains.txt"
endodata = processData( f1,skip=1,dataRowsLabels = seq( 16,20 ) )
endodata$endo_quant_log = endodata$data
endodata$endo_quant_log[ endodata$endo_quant_log<1 ] = endodata$endo_quant_log[ endodata$endo_quant_log<1 ] + 1
endodata$endo_quant_log = log2( endodata$endo_quant_log )
# rename columns
colnames( endodata )[ which( names( endodata ) == "data_label" ) ] = "time"
colnames( endodata )[ which( names( endodata ) == "data" ) ] = "endo_quant"

f2 = "./dynamic_metabolome_08052014/Exometabolome_ParentalStrains.txt"
exodata = processData(f2,skip=1,dataRowsLabels = seq(16,20) )
exodata$exo_quant_log = exodata$data
exodata$exo_quant_log[ exodata$exo_quant_log<1 ] = exodata$exo_quant_log[ exodata$exo_quant_log<1 ] + 1
exodata$exo_quant_log = log2( exodata$exo_quant_log )
# rename columns
names( exodata )[ which( names( exodata ) == "data_label" ) ] = "time"
names( exodata )[ which( names( exodata ) == "data" ) ] = "exo_quant"

f3 = "./dynamic_metabolome_08052014/growth_ParentalStrains.txt"
growthdata = processData(f3,skip=0,dataRowsLabels = seq(16,20) )
# rename attributes
levels(growthdata$metabolite)[1] <- "cellconc_1|ml"
levels(growthdata$metabolite)[2] <- "biovolume_ul|ml"
levels(growthdata$metabolite)[3] <- "singlecellvol_fl"
# rename columns
colnames( growthdata )[ which( names( growthdata ) == "metabolite" ) ] = "variable"
colnames( growthdata )[ which( names( growthdata ) == "data_label" ) ] = "time"
colnames( growthdata )[ which( names( growthdata ) == "data" ) ] = "value"
# cast to wide format
growthdata <- dcast( growthdata, strain + parent + time ~ variable, value.var = "value" )

# merge data
data = merge( endodata, exodata, all = TRUE )
data = merge( data, growthdata , all = TRUE )

# normalize using mixed model
# convert endometabolome data to matrix, use log2 data

# inititalize matrix
metabolites = unique(data$metabolite)
samples = sort(unique(paste(data$strain,data$time,sep="_")))
endo_matrix = matrix(0,nrow=length(metabolites),ncol=length(samples),dimnames=list(metabolites,samples))
exo_matrix = matrix(0,nrow=length(metabolites),ncol=length(samples),dimnames=list(metabolites,samples))

# fill matrices
for (i in 1:dim(data)[1]) {
	s = paste(data[i,]$strain,data[i,]$time,sep="_")
	m = levels(data[i,]$metabolite)[data[i,]$metabolite]
	v1 = data[i,]$endo_quant_log
	v2 = data[i,]$exo_quant_log
	endo_matrix[m,s] = v1
	exo_matrix[m,s] = v2
}
# assuming that each batch represents a bio rep
batch = factor(unlist(lapply(strsplit(colnames(endo_matrix),"_"),"[",1)))
cohort = factor(sapply(unlist(lapply(strsplit(colnames(endo_matrix),"_"),"[",1)),strtrim,1))

# remove metabolites with NAs
endo_matrix = endo_matrix[-which(apply(apply(endo_matrix,2,is.na),1,sum)>0),]
exo_matrix = exo_matrix[-which(apply(apply(exo_matrix,2,is.na),1,sum)>0),]

normalized_endo = normalizeMixed(endo_matrix,cohort=cohort,batch=batch)
normalized_exo = normalizeMixed(exo_matrix,cohort=cohort,batch=batch)

# add normalized data back
data = data %>%
	rowwise() %>%
	do({
		metabolite = levels(.$metabolite)[.$metabolite]
		strain = paste(.$strain,.$time,sep="_")
		if ( metabolite%in%rownames(normalized_endo) & strain%in%colnames(normalized_endo) ) {
			data.frame(.,endo_quant_log_normalized = normalized_endo[metabolite,strain])
		} else {
			data.frame(.,endo_quant_log_normalized = NA)
		}
	}) %>%
	do({
		metabolite = levels(.$metabolite)[.$metabolite]
		strain = paste(.$strain,.$time,sep="_")
		if ( metabolite%in%rownames(normalized_exo) & strain%in%colnames(normalized_exo) ) {
			data.frame(.,exo_quant_log_normalized = normalized_exo[metabolite,strain])
		} else {
			data.frame(.,exo_quant_log_normalized = NA)
		}
	})

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



data_complete_long = melt( data_complete )

# Save data and env ---------------------------------------------------

save( data, data_complete, data_complete_long, file="./dynamic_metabolome_08052014/data.RData" )

# Plot overall data trends ---------------------------------------------------

pdf("./dynamic_metabolome_08052014/metabolome.pdf", width=11.5,height=8)
# box plot of each metabolite
# endoboxplot <- ggplot(data = endodata, aes(x = metabolite, y = log(data))) +
# 	geom_violin(scale="width",fill="gray") +
# 	scale_x_discrete(name="") +
# 	scale_y_continuous(name=expression(paste("Level [ log10(",mu,"M) ]"))) +
# 	theme(axis.text.x=element_text(angle = -45, hjust = 0))

# Log2 metaboblites, Endo and Exo
ggplot(data = filter( data_complete_long, variable == "endo_quant_log_normalized")%>%filter(! is.na(value)), aes(x = time, y = value, group = parent,colour=parent, fill=parent)) +
  	geom_point() +
	geom_smooth() +
	facet_wrap( ~ metabolite, scales="free_y") +
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name=expression(paste("Level log2(",mu,"M)"))) +
    ggtitle("Endo- Metabolome Levels Across Strains") +
    scale_alpha_manual(values=c(.5,.25))
# Log2 metaboblites, Exo
ggplot(data = filter( data_complete_long, variable == "exo_quant_log_normalized" )%>%filter(! is.na(value)), aes(x = time, y = value, group = parent,colour=parent, fill=parent)) +
  	geom_point() +
	geom_smooth() +
	facet_wrap( ~ metabolite, scales="free_y") +
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name=expression(paste("Level log2(",mu,"M)"))) +
    ggtitle("Exo- Metabolome Levels Across Strains") +
    scale_alpha_manual(values=c(.5,.25))
# Rel t0 metaboblites, Endo
ggplot(data = filter( data_complete_long, variable == "endo_quant_rel" )%>%filter(! is.na(value)), aes(x = time, y = value, group = parent,colour=parent, fill=parent)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ metabolite, scales="free_y") +
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name="Relative Level") +
    ggtitle("Endo- Metabolome Levels Across Strains")
# Rel t0 metaboblites, Exo
ggplot(data = filter( data_complete_long, variable == "exo_quant_rel" )%>%filter(! is.na(value)), aes(x = time, y = value, group = parent,colour=parent, fill=parent)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ metabolite, scales="free_y") +
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name="Relative Level") +
    ggtitle("Exo- Metabolome Levels Across Strains")
# Rates metaboblites, Endo
d = filter( data_complete_long, variable == "endo_rate" )
d = d[ complete.cases(d), ]
ggplot(data = d, aes(x = time, y = value, group = parent,colour=parent, fill=parent)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ metabolite, scales="free_y") +
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name="Rate of Change (d[uM]/dt)") +
    ggtitle("Endo- Metabolome Rate of Change Across Strains")
# Rates metaboblites, Exo
d = filter( data_complete_long, variable == "exo_rate" )
d = d[ complete.cases(d), ]
ggplot(data = d, aes(x = time, y = value, group = parent,colour=parent, fill=parent)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ metabolite, scales="free_y") +
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name="Rate of Change (d[uM[/dt)") +
    ggtitle("Exo- Metabolome Rate of Change Across Strains")
# Rates metaboblites, Exo
ggplot(data = filter( data_complete_long, variable == "cellconc_1.ml" | variable == "biovolume_ul.ml" | variable == "singlecellvol_fl" ) %>% distinct(strain,time,variable), aes(x = time, y = value, group = parent, colour=parent)) +
	geom_point() +
	geom_smooth() +
	facet_wrap( ~ variable, scales="free_y") +
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name="Relative measurement")
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
