#' Cluster data
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return Clustering object of input data
#' @examples
#' cluster()
#' @export
#'
cluster = function(data, cluster_method = c("fuzzy")[1],
                   distance = c("cosine","euclidean")[1] ) {
  if (distance=="cosine") {
    dist_matrix = cosineDist(data)
  } else {
    dist_matrix = dist( data, method = distance )
  }
  # cluster
  if ( cluster_method == "fuzzy" ) {
    library(cluster)
    clustering = fanny( x = dist_matrix, k = 2 )
  }
  return(clustering)
}

#' Score clustering
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return Score representing distance between the actual clustering
#'  and the one proposed by model
#' @examples
#' score()
#' @export
#'
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

#' Permute
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return permuted scores
#' @examples
#' permute()
#' @export
#'
permute = function(genotype_vector,clusters,nresample = 10000,...) {
  # order clusters and genotypes_vector
  ecdf(unlist(mclapply(1:nresample,function(i){
    # can not do this globally as before because of biases
    # in genotype dist at some locations
    #genotype_vector = genotypes[sample(seq(1,dim(genotypes)[1]),1),]
    names(genotype_vector) = sample(names(genotype_vector))
    score(clusters,genotype_vector)
  })))
}

#' Cluster and Score
#'
#' The function \code{\link{cluster}}
#'
#' @param data
#' @return to_r
#' @examples
#' clustANDscore()
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
  data = data[which(apply(data,1,sum)>0,useNames=T),]
  # only use data for clustering
  data = data[intersect(colnames(genotypes),rownames(data)),]
  #covariate_clust = cluster(covariates)
  clustering = cluster(data,...)
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
    one_score = score(clustering$clustering,genotype_vector)
    permuted = permute(genotype_vector,clustering$clustering)
    matrix(c(one_score,permuted(one_score)),nrow=1,ncol=2,dimnames=list(rownames(genotypes)[x],c("score","pval")))
  }))
}
