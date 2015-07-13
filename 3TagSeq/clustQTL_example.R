library( GenomicRanges )
library( parallel )
library( ggbio )
options("mc.cores"=20)
source("/g/steinmetz/brooks/git/R-tools/multiplot.R")
# files 
load("/g/steinmetz/project/GenPhen/data/3tagseq/all/mergedCounts.rda")
load("/g/steinmetz/project/GenPhen/data/3tagseq/counts.rda")

load( "/g/steinmetz/brooks/genphen/qtl_endometabolome_23042015/geno_mrk.RData" )

g = makeGRangesFromDataFrame(tx_3utr,
        seqnames.field=c("seqnames"),
        start.field=c("start"),
        end.field=c("end"))

load("/g/steinmetz/brooks/3prime_Tfill/clust_qtls.rda")

# gene with sig qtl
p_matrix = do.call(cbind,lapply(clust_qtls,function(i){-log10(i[,2]+1e-4)}))
sig_qtls = names(which(apply(p_matrix,2,function(i)sum(i>3))>0))

#i = sample(sig_qtls,1)
i = "YIL014C-A"
# phenotype index
i_ind = which(tx_3utr$Name==i)
# 7024
data = t(as.matrix(mcols(subsetByOverlaps(cts,g[i_ind]))))
# need to clean up data names
rnames = sapply( sub("X","",sapply( rownames(data), function(i) strsplit(i,"_")[[1]][1] ) ), function(i) {
		if (nchar(i)==2){ 
			i = paste("0",i,sep="")
		} 
		i
	} )
rownames(data) = rnames
# remove insig peaks
peak_cutoff = sig_peaks( data )
data[ , colSums(data)<peak_cutoff ] = 0
# otherwise will bias clustering
# data = data[which(apply(data,1,sum)>0,useNames=T),]
# only use data with genotypes for clustering
data = data[intersect(colnames(geno),rownames(data)),]
clust_qtl = clustANDscore(data,geno,distance="binary")
qtls = clust_qtl
qtl_name = i
plotManhattan(clust_qtl,mrk, i,trx_annot=tx_3utr)

peaks = findQTLPeaks(clust_qtl,mrk)

peak_genotypes = do.call(rbind,lapply(names(peaks),function(i){getGenotypes(i,geno)}))
rownames(peak_genotypes) = names(peaks)

# test whether QTL looks "real"

# peak to analyze
x = 1

A = names(which(peak_genotypes[x,] == 1))
A = intersect( A, rownames(data) )
B = names(which(peak_genotypes[x,] == 2))
B = intersect( B, rownames(data) )


# correlation
all_cor = cor(t(data[c(A,B),]),use="pair",method="pearson")
in_a = all_cor[A,A]
in_b = all_cor[B,B]
out_ab = all_cor[A,B]
out_ba = all_cor[B,A]

#mean(in_a[upper.tri(in_a)],na.rm=T)
#mean(in_b[upper.tri(in_b)],na.rm=T)
#mean(out_ab[upper.tri(out_ab)],na.rm=T)
#mean(out_ba[upper.tri(out_ba)],na.rm=T)

boxplot(list(in_a[upper.tri(in_a)],in_b[upper.tri(in_b)],out_ab[upper.tri(out_ab)]))

# distance
all_dist = as.matrix(dist((data[c(A,B),])))
in_a = all_dist[A,A]
in_b = all_dist[B,B]
out_ab = all_dist[A,B]
out_ba = all_dist[B,A]

# mean(in_a[upper.tri(in_a)],na.rm=T)
# mean(in_b[upper.tri(in_b)],na.rm=T)
# mean(out_ab[upper.tri(out_ab)],na.rm=T)
# mean(out_ba[upper.tri(out_ba)],na.rm=T)

boxplot(list(in_a[upper.tri(in_a)],in_b[upper.tri(in_b)],out_ab[upper.tri(out_ab)]))


# heatmap
# only keep data with at least one positive value

clustering = cluster(data)
cluster1 = intersect(names(which(clustering$clustering==1)),colnames(geno))
cluster2 = intersect(names(which(clustering$clustering==2)),colnames(geno))

library(pheatmap)
pheatmap(data[c(A,B),],cluster_cols=F,cluster_rows=T,scale="row",annotation_row=data.frame(rbind(cbind(A,"S"),cbind(B,"Y")),row.names=1))


boxplot(list(log2(apply(data[A,],1,sum)),log2(apply(data[B,],1,sum))))
# sum(apply(p_matrix,2,function(i)sum(i>3))>0)
# 2823

# analyze expression component of 3 isoform clustering
# see if expression can explain all/most clustering
pb <- txtProgressBar(min = 0, max = dim(tx_3utr)[1], style = 3)
expression_covary = mclapply(seq(1,dim(tx_3utr)[1]),function(i_ind){
	setTxtProgressBar(pb, i_ind)
	#print(i_ind)
	to_r = try({
		data = t(as.matrix(mcols(subsetByOverlaps(cts,g[i_ind]))))
		# need to clean up data names
		rnames = sapply( sub("X","",sapply( rownames(data), function(i) strsplit(i,"_")[[1]][1] ) ), function(i) {
				if (nchar(i)==2){ 
					i = paste("0",i,sep="")
				} 
				i
			} )
		rownames(data) = rnames
		# remove insig peaks
		peak_cutoff = sig_peaks( data )
		data[ rowSums(data)>=peak_cutoff ,] = 0
		# otherwise will bias clustering
		data = data[which(apply(data,1,sum)>0,useNames=T),]
		# only use data for clustering
		data = data[intersect(colnames(geno),rownames(data)),]
		expression = rowSums(data)
		score(cluster(data)$clustering,cluster(expression)$clustering)
	},silent=T)
	if (class(to_r)=="try-error") {
		return(NULL)
	} else{
		return(to_r)
	}
})
close(pb)
names(expression_covary) = tx_3utr$Name

expression_linked_QTLs = sapply(sig_qtls,function(i){expression_covary[[i]]<min(clust_qtls[[i]][,1])})
# 93.8% of 3' clustQTLs are expression linked
sum(unlist(expression_linked_QTLs))/length(expression_linked_QTLs)
