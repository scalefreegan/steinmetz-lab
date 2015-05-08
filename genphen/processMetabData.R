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

# NOTE: this function is specific to data format supplied by Nicole, which looks like this:
# Intracellulare concentration [µM]												
# 	ISOCIT											AKG
# 	S1	S2	S3	S4	S5	S6	S7	Y1	Y2	Y3	Y4	S1
# 1	3844.50	3584.21	3947.04	4127.05	4275.45	4423.63	3901.31	2153.74	1961.84	2237.32	2114.77	34.1
# 2	2575.07	3084.54	3259.4	3199.84	3026.13	3414.68	3010.6	2151.62	2310.91	2146.4	1957.09	49.3
# 3	2579.87	2098.17	2552.03	3352.06	2661.11	2612.89	2529.99	1581.68	1648.46	1807.51	1589.14	31.8
# 4	2153.68	2157.96	2081.88	2047.92	2089.03	2078.23	1968.13	1341.41	1277.06	1430.45	1237.51	79.68
# 5	1609.98	1692.87	1625.89	1511.57	1674.03	1644.89	1489.67	1198.9	1207.58	1238.84	1175.75	29.49

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

f1 = "./dynamic_metabolome_08052014/Endometabolome_ParentalStrains.txt"
endodata = processData(f1,skip=1,dataRowsLabels = seq(16,20) )
pdf("./dynamic_metabolome_08052014/endometabolome.pdf", width=11.5,height=8)
ggplot(data = endodata, aes(x = data_label, y = data, group = parent,colour=parent)) + 
	geom_point() + 
	geom_smooth() + 
	facet_wrap( ~ metabolite, scales="free_y") + 
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name=expression(paste("Level (",mu,"M)")))
dev.off()

f2 = "./dynamic_metabolome_08052014/Exometabolome_ParentalStrains.txt"
exodata = processData(f2,skip=1,dataRowsLabels = seq(16,20) )
pdf("./dynamic_metabolome_08052014/exometabolome.pdf", width=11.5,height=8)
ggplot(data = exodata, aes(x = data_label, y = data, group = parent,colour=parent)) + 
	geom_point() + 
	geom_smooth() + 
	facet_wrap( ~ metabolite, scales="free_y") + 
	scale_x_discrete(name="Time (hr)") +
    scale_y_continuous(name=expression(paste("Level (",mu,"M)")))
dev.off()

