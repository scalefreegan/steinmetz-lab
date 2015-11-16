library(shiny)
library(clustQTL)
library(dplyr)
library(reshape2)
library(DT)
library(parallel)
library(GenomicRanges)
library(ggplot2)
library(funqtl)
library(stringr)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("org.Sc.sgd.db")

library(RCurl)
library(httr)
set_config( config( ssl_verifypeer = 0L ) )
# Global variables ---------------------------------------------------

id2name = id2name(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
type = "mlod"

# gene name map
x = org.Sc.sgdGENENAME
keys <- mappedkeys(x)
# Convert to a list
gname <- as.list(x[keys])

# short description map
x = org.Sc.sgdALIAS
keys <- mappedkeys(x)
# Convert to a list
dname <- as.list(x[keys])

# short description map
x = org.Sc.sgdDESCRIPTION
keys <- mappedkeys(x)
# Convert to a list
dname_long <- as.list(x[keys])

# Web resources ---------------------------------------------------
library(funqtl)
library(GenomicRanges)
library(rtracklayer)
library(jsonlite)

DDIR = "/g/steinmetz/brooks/genphen/metabolome/qtls"   

load(file.path(DDIR,"endometabolite_full_12102015.rda"))
load(file.path(DDIR,"genotypes_S288c_R64.rda"))
load(file.path(DDIR,"mQTLs_comball_funqtl_2014.rda"))

data = lapply(seq(1,length(mQTLs_funqtl_2014)), function(i){
        if (class(mQTLs_funqtl_2014[[i]])=="try-error") {
            return(NULL)
        } else {
            #print(i)
            o = list()
            o[["qtl"]] = mQTLs_funqtl_2014[[i]][["qtls_alt"]]
            o[["permout"]] = mQTLs_funqtl_2014[[i]]$permout
            return(o)
        }
    })
names(data) = names(mQTLs_funqtl_2014)

# remove NULLs
data = data[!sapply(data,is.null)]

################################################################################
## 
outfolder = "/g/steinmetz/project/GenPhen/web/mQTL_plots/20151113"
if (!file.exists(outfolder))  dir.create(outfolder)

fdrRange = seq(0, 50, by=1)
traits = names(data)

plot_tab = data.frame(
    FDR = rep(fdrRange, times = length(traits)),
    metabolite = rep(traits, each = length(fdrRange)), 
    fname = as.character(NA),
    stringsAsFactors=FALSE
)

library(Cairo)
for(thisTrait in traits) {
    
    cat(thisTrait, "\n")
    qtls = data[[thisTrait]]$qtl
    names(qtls)[names(qtls) == type] = "pval"
    qtls[,"pval"]  = 10^-qtls[,"pval"]    
    
    alpha_5 = summary(data[[thisTrait]]$permout[,type],.05)
    subfolder = file.path(outfolder, thisTrait)
    if (!file.exists(subfolder))  dir.create(subfolder)
    
    for(thisFDR in fdrRange){
        alpha_10 = summary(data[[thisTrait]]$permout[,type],thisFDR/100)
        
        fname =  file.path(subfolder, paste0("mQTL_", thisTrait, "_FDR_", thisFDR,".png"))        
        ymax = max(max(alpha_5, alpha_10), max(-log10(qtls[,"pval"]))) 
        
        CairoPNG(fname, width=1500, height=600)
        clustQTL::plotManhattan(qtls, mrk, qqman = TRUE, show = TRUE,
            suggestiveline = alpha_5, genomewideline = alpha_10, 
            ylab = "LOD", ylim = c(0,ymax+2))
        legend("topright", y.leg[i], c("5% FDR",paste(rownames(alpha_10)[1],"FDR")), 
            lty = c(1, 1), col = c("blue", "red"))
        dev.off()
        plot_tab[plot_tab$FDR == thisFDR & plot_tab$metabolite == thisTrait,"fname"] = fname    
    }
} 

save(plot_tab, file=file.path(DDIR, "plot_tab.rda"))
