library( GenomicRanges )
library( parallel )
options("mc.cores"=20)
# files
load("/g/steinmetz/project/GenPhen/data/3tagseq/all/mergedCounts.rda")
load("/g/steinmetz/project/GenPhen/data/3tagseq/counts.rda")

load( "/g/steinmetz/brooks/genphen/qtl_endometabolome_23042015/geno_mrk.RData" )


normalize = function(x) {
	t(apply(x,1,function(i)(i-mean(i))/sd(i)))
}

g = makeGRangesFromDataFrame(tx_3utr,
        seqnames.field=c("seqnames"),
        start.field=c("start"),
        end.field=c("end"))

sig_peaks = function(peakMatrix, nresamples = 1000, alpha = 0.01) {
	# resample peak matrix
	# matrix should be nxm matrix where
	# n = position
	# m = observation (eg strain)
	# each value at pos n,m is count
	peakSamples = unlist(lapply(1:nresamples,function(i){
		rowSums(do.call(cbind,lapply(1:dim(peakMatrix)[2],function(i){sample(peakMatrix[,i])})))
	}))
	# return cutoff peak height at threshold
	quantile(peakSamples,1-alpha)
}

# rep

# strainwise isoform correlation, replicates removed
n = dim(mcols(g))[1]
strain_isoform_corr = unlist(mclapply(1:n,function(i){
	if (i%%1000==0){
		cat(paste(round((i/n)*100),"%\n",sep=""))
	}
	g_cts = as.matrix( mcols( subsetByOverlaps(cts,g[i]) ) )
	# filter peaks
	cts_t = sig_peaks( g_cts )
	# set counts at insig peaks to zero
	g_cts[ rowSums(g_cts)<cts_t ,] = 0
	gs_cor = try(cor(g_cts),silent=T)
	# make replicate mask
	rstrain = sapply(rownames(gs_cor),function(i){strsplit(i,split="_")[[1]][1]})
	cstrain = sapply(colnames(gs_cor),function(i){strsplit(i,split="_")[[1]][1]})
	repmask = do.call(rbind,lapply(rstrain,function(x){cstrain%in%x}))
	gs_cor[repmask] <- NA
	gs_cor = gs_cor[upper.tri(gs_cor)]
	# remove nas
	gs_cor = gs_cor[!is.na(gs_cor)]
	return(gs_cor[!is.na(gs_cor)])
	}))

# generate evernote plots
pdf("/g/steinmetz/brooks/3prime_Tfill/peakfiltering.pdf")
for (i in 1:100) {
	print(i)
	n = tx_3utr[i,"Name"]
	g_cts = as.matrix( mcols( subsetByOverlaps(cts,g[i]) ) )
	cts_t = sig_peaks( g_cts )
	g_cts2 = g_cts
	g_cts2[ rowSums(g_cts)<cts_t ,] = 0
	try({par(mfrow=c(2,1))
	matplot((g_cts),type="l",xlab="Position relative to TTS",ylab="Count",main=paste(n,"3' isoforms before filtering",sep=": ") )
	matplot((g_cts2),type="l",xlab="Position relative to TTS",ylab="Count",main=paste(n,"3' isoforms after filtering",sep=": ") )})
}
dev.off()

pdf("/g/steinmetz/brooks/3prime_Tfill/strain_isoform_corr.pdf")
	hist(strain_isoform_corr,breaks=1000,main="3' isoform similarity across strains, replicates removed",xlab="Pearson correlation")
dev.off()

# clustQTL



namematch = function(string, stringVector) {
	library(stringdist)
	match_ind = amatch(string, stringVector,maxDist=nchar(string)-mean(sapply(stringVector,nchar))+1)
	if (!is.na(match_ind)) {
		return( stringVector[ match_ind ] )
	} else {
		return( match_ind )
	}
}

cluster = function(data, cluster_method = c("fuzzy")[1], distance_metric = c("euclidean")[1] ) {
	dist_matrix = dist( data, method = distance_metric )
	# cluster
	if ( cluster_method == "fuzzy" ) {
		library(cluster)
		clustering = fanny( x = dist_matrix, k = 2 )
	}
	return(clustering)
}

score = function(clusters, genotype_vector,score_method = c("hamming")[1], verbose=FALSE) {
	# order clusters and genotypes_vector
	clusters = clusters[intersect(names(genotype_vector),names(clusters))]
	genotype_vector = genotype_vector[names(clusters)]
	# find best match to assign genotypes
	genotype_vector = as.numeric(as.factor(genotype_vector))
	if ( cor(genotype_vector,clusters) < 0 ) {
		# switch labels
		new_genotype_vector = genotype_vector
		new_genotype_vector[genotype_vector==1] = 2
		new_genotype_vector[genotype_vector==2] = 1
		genotype_vector = new_genotype_vector
	}
	# compute hamming distance
	if (score_method == "hamming") {
		if (verbose){
			cat("Scoring by Hamming distance\n")
		}
		score = sum(clusters!=genotype_vector)
	}
	return(score)
}

permute = function(genotypes,clusters,nresample = 10000,...) {
	# order clusters and genotypes_vector
	ecdf(unlist(mclapply(1:nresample,function(i){
		genotype_vector = genotypes[sample(seq(1,dim(genotypes)[1]),1),]
		names(genotype_vector) = sample(names(genotype_vector))
		score(clusters,genotype_vector)
	})))
}

clustANDscore = function(data, genotypes,...) {
	# change data names to match genotypes
	# new_names = sapply(rownames(data),namematch,colnames(genotypes))
	# remove names that don't match from data
	# data = data[!is.na(new_names),]
	# rownames(data) = new_names[!is.na(new_names)]
	# only keep data with at least one positive value
	# otherwise will bias clustering
	data = data[which(apply(data,1,sum)>0,useNames=T),]
	# only use data for clustering
	data = data[intersect(colnames(genotypes),rownames(data)),]
	#covariate_clust = cluster(covariates)
	clustering = cluster(data)
	permuted = permute(genotypes,clustering$clustering)
	genotype_template = rownames(data)
	names(genotype_template) = genotype_template
	to_r = do.call(rbind,mclapply(seq(1:dim(genotypes)[1]),function(x){
		genotype_vector=sapply(genotype_template,function(i){
			if (class(try(genotypes[x,i],silent=T))!="try-error"){
				return(genotypes[x,i])
			} else {
				return(NA)
			}
		})
		genotype_vector = genotype_vector[!is.na(genotype_vector)]
		one_score = score(clustering$clustering,genotype_vector)
		matrix(c(one_score,permuted(one_score)),nrow=1,ncol=2,dimnames=list(rownames(genotypes)[x],c("score","pval")))
	}))
}

format4manhattan = function( qtls, mrk, nresample = 10000 ) {
	mrk$p = -log10( qtls[,"pval"]+ 1/nresample )
	return(mrk)
}

plotManhattan = function( qtls, mrk, qtl_name = "", trx_annot = NULL, cutoff = 3, gene_annot_range = c(5000,5000) ) {
	library(ggbio)
	mrk2 = format4manhattan( qtls,mrk )
	if ( (qtl_name!="") & ( !is.null(trx_annot) ) ) {
		# annotate gene
		trx_info = trx_annot[ which(trx_annot$Name == qtl_name), ]
		trx_granges = GRanges(seqnames=trx_info$seqnames,ranges=IRanges(trx_info$start-gene_annot_range[1],trx_info$end+gene_annot_range[2]))
		names(trx_granges) = qtl_name
		plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,ylab="-log10(pval)",main=qtl_name,ylim=c(0,4.5),highlight.gr = trx_granges)
	} else {
		plotGrandLinear(mrk2, aes(y = p),spaceline = TRUE,cutoff=cutoff,ylab="-log10(pval)",main=qtl_name,ylim=c(0,4.5))
	}

	#return(mrk2)
}

findQTLPeaks = function( qtls, mrk, cutoff = 3 ) {
	# find top two peaks
	mrk2 = format4manhattan( qtls, mrk )
	candidates = which( mrk2$p>cutoff)
	ind = 1
	final_candidates = cbind(candidates[ind],mrk2$p[candidates[ind]])
	ind = ind + 1
	while (ind<=length(candidates)) {
		if (diff(candidates)[ind-1]<10) {
			# single peak
			# only replace val if higher
			if ( mrk2$p[candidates[ind]] > final_candidates[ dim(final_candidates)[1],2] ) {
				final_candidates[ dim(final_candidates)[1],] = cbind(candidates[ind],mrk2$p[candidates[ind]])
			}
		} else {
			final_candidates = rbind(final_candidates,cbind(candidates[ind],mrk2$p[candidates[ind]]))
		}
		ind = ind + 1
	}
	to_r = mrk2[paste("mrk",final_candidates[,1],sep="_"),]
	return(to_r)
}

getGenotypes = function(marker,geno) {
	geno[marker,]
}

# try it
pb <- txtProgressBar(min = 0, max = dim(tx_3utr)[1], style = 3)
clust_qtls = lapply(seq(1,dim(tx_3utr)[1]), function(i){
#clust_qtls = lapply(seq(5,6), function(i){
	setTxtProgressBar(pb, i_ind)
	data = t(as.matrix(mcols(subsetByOverlaps(cts,g[i]))))
	# need to clean up data names
	rnames = sapply( sub("X","",sapply( rownames(data), function(i) strsplit(i,"_")[[1]][1] ) ), function(i) {
			if (nchar(i)==2){
				i = paste("0",i,sep="")
			}
			i
		} )
	rownames(data) = rnames
	try(clustANDscore( data, geno ) )
	})
close(pb)
names(clust_qtls) = tx_3utr$Name
# clean
clust_qtls = clust_qtls[!sapply(clust_qtls,class)=="try-error"]
save(clust_qtls,file="/g/steinmetz/brooks/3prime_Tfill/clust_qtls.rda")

load("/g/steinmetz/brooks/3prime_Tfill/clust_qtls.rda")
