#' Cluster data
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return Clustering object of input data
#' @export
#'
clustqtl_cluster = function(data, cluster_method = c("fuzzy")[1],
                   distance = c("cosine","euclidean")[1] ) {
  if (distance=="cosine") {
    dist_matrix = cosineDist(data)
  } else {
    dist_matrix = dist( data, method = distance )
  }
  # cluster
  if ( cluster_method == "fuzzy" ) {
    clustering = cluster::fanny( x = dist_matrix, k = 2 )
  }
  return(clustering)
}

#' Score clustering
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return Score representing distance between the actual clustering
#'  and the one proposed by model using binomial distribution
#' @export
#'
clustqtl_score = function(clusters, genotype_vector,verbose=FALSE) {
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
  # compute pval by exact binomial test
  # calc hypothesized prob success given
  # genotype vector and cluster composition
  # p(1,1) + p(2,2)
  p_hyp = sum(clusters==1)/length(clusters)*sum(genotype_vector==1)/length(genotype_vector) +
    sum(clusters==2)/length(clusters)*sum(genotype_vector==2)/length(genotype_vector)
  score = binom.test(sum(clusters==genotype_vector),length(clusters),p_hyp,"greater")
  return(score)
}

#' Cluster and Score
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return to_r
#' @export
#'
#'
clustANDscore = function(data, genotypes,...) {
  # change data names to match genotypes
  # new_names = sapply(rownames(data),namematch,colnames(genotypes))
  # remove names that don't match from data
  # data = data[!is.na(new_names),]
  # rownames(data) = new_names[!is.na(new_names)]
  # only keep data with at least one positive value
  # otherwise will bias clustering
  # data = data[which(apply(data,1,sum)>0,useNames=T),]
  # only use data for clustering
  data = data[intersect(colnames(genotypes),rownames(data)),]
  #covariate_clust = cluster(covariates)
  clustering = clustqtl_cluster(data,...)
  #permuted = permute(genotypes,clustering$clustering)
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
    one_score = clustqtl_score(clustering$clustering,genotype_vector)
    o = cbind(one_score$statistic[[1]],one_score$p.value)
    colnames(o) = c("matches","pval")
    return(o)
  }))
  to_r[,"pval"] = p.adjust(to_r[,"pval"],method = "BH")
  return(to_r)
}
