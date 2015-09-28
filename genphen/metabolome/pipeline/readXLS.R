#! /usr/bin/env Rscript
#
#-------------------------------------------------------------------#
# Process Metabolome Data from XLS spreadsheet
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

library(xlsx)

# Expected data format ---------------------------------------------------

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

processData = function( f1, f2, startRow = 2, mLoc = "endo", isParental = "false", ... ) {
	# f1 is relative time - defined by "cultivation phase"
	# eg. Endometabolome_1B_25A_sorted by cultivation phase.xlsx
	# f2 contains measurements at absolute time
	# eg. Endometabolome_25A_sorted by cultivation time.xlsx
	wb1 = loadWorkbook(f1)
	wb2 = loadWorkbook(f2)
	sheets = intersect(names(getSheets(wb1)),names(getSheets(wb2)))
	pb <- txtProgressBar(min = 0, max = length(sheets), style = 3)
	o = do.call(rbind,lapply((1:length(sheets)),function(i){
		#print(i)
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
					fpath = f1
				} else {
					touse = thisf2
					fpath = f2
				}
				time = touse[,1]
				time_format = k
				value = abs(as.numeric(touse[,j]))
				# if (is.na(value[1])) {
				#   time = time[2:length(time)]
				#   value = value[2:length(value)]
				# }
				#value.log2 = value
				#value.log2[value>0 & !is.na(value)] = log2(value.log2)
				#relative.log2 = sapply(value.log2,function(i){i/value.log2[1]})
				#derivative.log2 = c(NA,diff(value.log2)/diff(time))
				data.frame(
					fpath,
					metabolite,
					mLoc,
					strain,
					replicate,
					isParental,
					time,
					time_format,
					value)
			}))
		}))
		# remove NA
		# to_r = to_r[!is.na(to_r$value),]
		return(to_r)
	}))
	close(pb)
  # return a R data.frame
	return(o)
}
