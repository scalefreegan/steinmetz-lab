library( GenomicRanges )
library( parallel )
options("mc.cores"=20)
# files 
load("/g/steinmetz/project/GenPhen/data/3tagseq/all/mergedCounts.rda")
load("/g/steinmetz/project/GenPhen/data/3tagseq/counts.rda")

normalize = function(x) {
	t(apply(x,1,function(i)(i-mean(i))/sd(i)))
}

g = makeGRangesFromDataFrame(tx_3utr,
        seqnames.field=c("seqnames"),
        start.field=c("start"),
        end.field=c("end"))

n = dim(mcols(g))[1]
cluster_btwnss_2mer = unlist(mclapply(1:n,function(i){
	if (i%%1000==0){
		cat(paste(round((i/n)*100),"%\n",sep=""))
	}
	tmp = subsetByOverlaps(cts,g[i])
	km = try(kmeans(t(as.matrix(mcols(tmp))),2),silent=T)
	if (class(km)!="try-error"){
		return(abs(diff(km$size)))
		#km$betweenss
	} else {
		return(NULL)
	}
}))


pheatmap(t(as.matrix(mcols(subsetByOverlaps(cts,g[i])))),clust_cols=F)
matplot(rowSums(as.matrix(mcols(subsetByOverlaps(cts,g[6])))),type="l")

m = as.matrix(mcols(subsetByOverlaps(cts,g[6])))
tmp = lapply(1:1000,function(j){
	rowSums(do.call(cbind,lapply(1:dim(m)[2],function(i){sample(m[,i])})))
})
matplot(tmp,type="l")



tmp = cor(as.matrix(mcols(subsetByOverlaps(cts,g[7]))))
hist(tmp[upper.tri(tmp)])

# strainwise isoform correlation
n = dim(mcols(g))[1]
n = 100
stain_isoform_corr = unlist(mclapply(1:n,function(i){
	if (i%%1000==0){
		cat(paste(round((i/n)*100),"%\n",sep=""))
	}
	g_cts = subsetByOverlaps(cts,g[i])
	gs_cor = try(cor(as.matrix(mcols(g_cts))),silent=T)
	gs_cor = gs_cor[upper.tri(gs_cor)]
	}))