# after getting the poly A site counts using 
# bedtool with genomeCoverageBed
# combine everything into a count table
# Author: czhu
###############################################################################
## read the counts tables
myfunc = function(thisfolder) {
    require(rtracklayer)
    ## use unionBedGraphs to combine the bed files
    ## !!!the output from bedtools doesn't have the # in the header
    ## you need to add it manually!!!
    rv = lapply(c("plus.bedGraph","minus.bedGraph"), function(filename){
            infile = file.path(thisfolder, filename)
            strnd = ifelse(grepl("plus", filename),"+","-")
            bed = import(infile, "bedGraph", trackLine=TRUE)
            strand(bed) = strnd
            nm = basename(strsplit(readLines(infile,n=1), "\t")[[1]][-(1:3)])
            nm = sapply(strsplit(nm,"\\."),"[",1)
            #nm = gsub(paste0(ifelse(strnd=="+","_plus","_minus"),".bedGraph"),"",nm)        
            colnames(values(bed)) = nm
            bed                 
        })
    
    cts  = do.call(c,rv)
    cts    
}

folder = "/tmpdata/czhu/GenPhen/3tagseq/counts"
cts = myfunc(folder)
save(cts, file=file.path(folder, "rawcounts.rda"))

require("GenomicRanges")
nm = colnames(values(cts))
## extract metadata
## example:
# C6LRUACXX_SxY_Tagseq_P1B_15s009415-2-1_Tekkedil_lane843Abiorep1_sequence
pd = data.frame(
    run = sub("^(.+?)_.+Tagseq_(P..)_.+_(lane.)(.+?)(biorep.+?)_.+?$","\\1", nm),
    libraryPlate = sub("^(.+?)_.+Tagseq_(P..)_.+_(lane.)(.+?)(biorep.+?)_.+?$","\\2", nm),
    lane = sub("^(.+?)_.+Tagseq_(P..)_.+_(lane.)(.+?)(biorep.+?)_.+?$","\\3", nm),
    strain = sub("^(.+?)_.+Tagseq_(P..)_.+_(lane.)(.+?)(biorep.+?)_.+?$","\\4", nm),
    biorep = sub("^(.+?)_.+Tagseq_(P..)_.+_(lane.)(.+?)(biorep.+?)_.+?$","\\5", nm) 
    ,stringsAsFactors=FALSE)

pd$strain = sub("0$","",pd$strain)

## merge the technical replicates
f= pd$sample = paste( pd$strain, pd$biorep, sep="_")
mat = as.matrix(values(cts))

mmat = do.call(cbind,mclapply(unique(f), function(thisF){
            wh = which(f == thisF)
            if(length(wh)==1){
                return(mat[,wh])
            } else{
                return(rowSums(mat[,wh]))
            }
        },mc.cores=10))
colnames(mmat) = unique(f)

values(cts) = mmat
folder="."
save(cts,pd,file=file.path(folder,"mergedCounts.rda"))




csums = colSums(mat)
mat[mat==0] = NA

plotFolder = "/g/steinmetz/project/GenPhen/data/3tagseq/plots"
if (!file.exists(plotFolder))  dir.create(plotFolder)

library(Cairo)
#CairoPNG(file.path(plotFolder, "distr_overall_size.png"), width=1500, height=500)
CairoPDF(file.path(plotFolder, "raw_distr_overall_size.pdf"), width=10, height=8)
#barplot(log10(csums),las=2)
hist(log10(csums),100, main="Distribution of sample size", 
    xlab="Number of mappable reads")
dev.off()

CairoPDF(file.path(plotFolder, "raw_distr_yield_by_lane.pdf"), width=12, height=12)
par(mfrow=c(3,3))
f = paste(pd$run, pd$lane)
for(thisF in unique(f)){
    hist(log10(csums)[f==thisF],100, main=paste("Distribution of sample size", thisF), 
        xlab="Number of mappable reads")    
}
dev.off()

## technical replicates
library(LSD)
f = paste(pd$strain, pd$biorep, sep="_")
subPlotFolder = file.path(plotFolder, "raw_tech_rep")
if (!file.exists(subPlotFolder))  dir.create(subPlotFolder)
subf = names(table(f)[table(f)>=2])
for(thisF in subf){
    CairoPNG(file.path(subPlotFolder, paste0(thisF,".png")), width=1200, height=1200)
    submat = mat[,f==thisF]
    submat = log2(submat[rowSums(submat>=5)>0,])
    heatpairs(submat)    
    dev.off()
}


