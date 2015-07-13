library(MASS)
library( GenomicRanges )

set.seed(1)

# files 
# load("/g/steinmetz/project/GenPhen/data/3tagseq/all/mergedCounts.rda")
# load("/g/steinmetz/project/GenPhen/data/3tagseq/counts.rda")

# fit properties of transcriptome data

# mean expression level
y = rowMeans(normCounts)
fit_e <- fitdistr(y,"Poisson")
# 172.18 for transcriptome data
#    lambda  
#  172.17929 
# (  0.15387)

# length of alt 3' window
y = tx_3utr$width
fit_l <- fitdistr(y,"normal")
#     mean        sd   
#  366.2071   131.1248 
# (  1.5377) (  1.0873)

# number of alt ends per region
g = makeGRangesFromDataFrame(tx_3utr,
        seqnames.field=c("seqnames"),
        start.field=c("start"),
        end.field=c("end"))

tmp = subsetByOverlaps(g,cts)

# expression level
e = round( rpois(1,172.18) )
# phe length
l = round( rnorm(1,366.2,131.1) )
phe = rep(0,l)

# 
p = table(ceiling(runif( e, 0,1 )*l))


sapply( runif( e, 0,1 ), function(i){ rpois(1,i) })

genotype