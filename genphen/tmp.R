
runPCA <- function(mat = 'Unadjusted matrix') eigen(cov(apply(mat, 2, function(i) (i - mean(i))/sd(i))))
pca <- runPCA(metabolome_data_mixednorm)

afterPCA <- function(
 matAdjust = 'Centered matrix',
 meanList = 'List of column means of original (unadjusted) matrix',
 eigenList = 'List of eigenvalues and eigenvectors of adjust matrix covariance matrix',
 n = 'selected PC\'s to remove',
 specific_select = 'If True: n == 1:n, if False: just n\'th columns') {
 
 if (length(n) > ncol(matAdjust)) stop('N is higher than the number of PC\'s')
 if (!specific_select & length(n) > 1) stop('Use a single number when selecting up to n\'th PC')
 if (!specific_select) n <- 1:n
 
 to_r = t(eigenList$vectors[,-n] %*% (t(eigenList$vectors[,-n]) %*% t(matAdjust))) + (matrix(meanList, nrow = nrow(matAdjust), ncol = ncol(matAdjust),byrow=T))
 colnames(to_r) = colnames(matAdjust) 
 return(to_r)

}
 
reconstMatrix <- afterPCA(
 matAdjust = apply(metabolome_data_mixednorm, 2, function(i) i - mean(i)),
 meanList = apply(metabolome_data_mixednorm, 2, mean),
 eigenList = pca,
 n = 1:15,
 specific_select = TRUE
)
