# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#
# Author: Alexandra Jauhiainen, alexandra.jauhiainen@ki.se

updateMean <- function(X,y,Z,Sigma,samples,return.fit=F){
  # Internal function to update mean in the iterative normalization procedure.
  # Args:
  #   X,Z: regressor matrices
  #   y: response vector
  #   Sigma: current estimate of covariance matrix
  #   samples: vector of sample indicator
  #   return.fit: If TRUE the fitted mixed model is returned, in addition to the other parameters
  # Returns:
  #   Updated fitted values and updated residuals, and, optional, the fitted model.
  
  Sigma.inv <- solve(Sigma)
  
  n <- length(levels(samples))
  m <- length(y)/n
  
  e <- eigen(Sigma.inv)
  V <- e$vectors
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  Ytmp <- matrix(y,nrow=n,ncol=m)
  
  Y2 <-NULL
  for(i in 1:n){
    Y2 <- rbind(Y2,t(B%*%(Ytmp[i,])))
  }
  
  X2 <- matrix(NA,nrow=nrow(X),ncol=ncol(X))
  Z2 <- matrix(NA,nrow=nrow(Z),ncol=ncol(Z))
  for(samp in levels(samples)){
    pick <- samples==samp
    Xtmp <- X[pick,]
    Ztmp <- Z[pick,]
    X2[pick,]<-B%*%Xtmp
    Z2[pick,]<-B%*%Ztmp
  }
  # new values in Y2, X2, Z2
  fit3 <- hglm(X = X2, y = as.numeric(Y2), Z = Z2, family = gaussian(link = identity),conv=1e-6)
  fitted.vals <- matrix(fit3$fv,nrow=n,ncol=m)
  resid.vals <- matrix(as.numeric(Y2)-fit3$fv,nrow=n,ncol=m)
  
  fitted.new <-NULL
  for(i in 1:n){
    fitted.new <- rbind(fitted.new,t(solve(B)%*%(fitted.vals[i,])))
  }
  
  resid.new <-NULL
  for(i in 1:n){
    resid.new <- rbind(resid.new,t(solve(B)%*%(resid.vals[i,])))
  }
  if(!return.fit)
    list(fitted=fitted.new,resid=resid.new)
  else
    list(fitted=fitted.new,resid=resid.new,fit=fit3)
}

updateSigma <- function(fitted.vals,Y,resid.vals,use.resid=TRUE,eps=eps){
  # Internal function to update the covariance matrix in the iterative normalization procedure.
  # Args:
  #   fitted.vals: Current fitted values.
  #   Y: Response values, as a matrix.
  #   resid.vals: Current residuals.
  #   use.resid: If TRUE the current residuals are used to estimate covaraince matrix. If FALSE the residuals a re-calculated.
  #   eps: small value to avoid numerical instability in estimating the covaraince matrix.
  # Returns:
  #   An updated estimate of the covariance matrix.
  
  m <- ncol(resid.vals)
  n <- nrow(resid.vals)
  
  M <- matrix(0,ncol=m,nrow=m)
  if(use.resid){
    for(i in 1:n){
      M<- M + resid.vals[i,] %*% t(resid.vals[i,])
    }
  }else{
    for(i in 1:n){
      M<- M + (Y-fitted.vals)[i,] %*% t((Y-fitted.vals)[i,])
    }
  }
  Sigma.cur <- M/n + eps*diag(m)
  Sigma.cur
}

normalizeMixed <- function(Y,cohort,batch,iterative=T,return.cov=F,scale.bv=F,conv.eps=1e-5){
  # Performs a mixed model normalization with a simultaneous estimation of
  # a covariance matrix by using an iterative procedure.
  # Args:
  #   Y: Response matrix, metabolites x samples.
  #   cohort: A vector indicating cohorts for each sample, can be set to NULL if no cohorts are defined.
  #   batch: A vector indicating batches for each sample, can be set to NULL if no batches are defined.
  #   iterative: if iteration is to be performed, defaults to true
  #   return.cov: If TRUE the estimated covariance matrix is returned together with the normalized values, provided that iterative is set to TRUE
  #   scale.bv: If TRUE, a model is fit initially to correct for differing variances for the batches. Requires the nlme package. 
  #   conv.eps: Small value to indicate convergence.
  # Returns:
  #   The normalized Y matrix, and, optional, the covariance matrix.
  require(hglm)
  
  Y <- t(Y)
  m <- ncol(Y)
  n <- nrow(Y)
  
  met <- sprintf("%s%2$0*3$d", "m", 1:m,nchar(m))
  met <- rep(met,each=n)
  sample <- sprintf("%s%2$0*3$d", "s", 1:n,nchar(n))
  sample <- rep(sample,times=m)
  
  if(is.null(batch) & is.null(cohort)){
    message("A model with only samples is used.")
    dat <- data.frame(y=as.numeric(Y),met=as.factor(met),sample=as.factor(sample))
    X<-model.matrix(~1+met,data=dat)
    Z <- model.matrix(~0+sample,data=dat)
  }
  if(is.null(cohort) & !is.null(batch)){
    message("A model with batches and samples is used.")
    batch <- as.factor(batch)
    batch <- rep(batch,times=m)
    dat <- data.frame(y=as.numeric(Y),batch=as.factor(batch),met=as.factor(met),sample=as.factor(sample))
    X<-model.matrix(~1+met,data=dat)
    Z <- model.matrix(~0+batch,data=dat)
    Z <- cbind(model.matrix(~0+sample,data=dat),Z)
    if(scale.bv){
      require(nlme)
      fit1 <- lme(y~met, random=~1|batch/sample,weights = varIdent(form=~1|batch),data=dat)
      dat$y <- residuals(fit1,type="pearson") + fitted(fit1)
    }
  }
  if(!is.null(cohort) & !is.null(batch)){
    message("A model with cohorts, batches, and samples is used.")
    cohort <- as.factor(cohort)
    cohort <- rep(cohort, times=m)
    dat <- data.frame(y=as.numeric(Y),cohort=as.factor(cohort),batch=as.factor(batch),met=as.factor(met),sample=as.factor(sample))
    X<-model.matrix(~1+met,data=dat)
    Z <- model.matrix(~0+cohort,data=dat)
    Z <- cbind(model.matrix(~0+batch,data=dat),Z)
    Z <- cbind(model.matrix(~0+sample,data=dat),Z)
    if(scale.bv){
      require(nlme)
      fit1 <- lme(y~met, random=~1|cohort/batch/sample,weights = varIdent(form=~1|batch),data=dat)
      dat$y <- residuals(fit1,type="pearson") + fitted(fit1)
    }
  }
  
  y = dat$y
  fit2 <- hglm(X = X, y = y, Z = Z, family = gaussian(link = identity))
  
  if(iterative){
    fitted.cur <- fitted.start <- matrix(fit2$fv,nrow=n,ncol=m)
    eps <- 1e-10
    Sigma.cur <- Sigma.start<-cov(matrix(dat$y-fit2$fv,nrow=n,ncol=m)) + eps*diag(m)
    
    j<-1
    converged<-F
    while(!converged){
      
      out.mean <- updateMean(X,y,Z,Sigma.cur,dat$sample)
      Sigma.new <- updateSigma(out.mean$fitted,Y,out.mean$resid,use.resid=T,eps=eps)
      
      diff<-abs(Sigma.cur[upper.tri(Sigma.cur)]-Sigma.new[upper.tri(Sigma.new)])
      if(all(diff<conv.eps)){
        converged<-T
        final.fit <- updateMean(X,y,Z,Sigma.new,dat$sample,return.fit=T)
      }
      Sigma.cur<-Sigma.new
      colnames(Sigma.cur)<-rownames(Sigma.cur)<-colnames(Y)
      j<-j+1
    }
    
    colnames(Sigma.cur)<-rownames(Sigma.cur)<-colnames(Y)
    colnames(Sigma.start)<-rownames(Sigma.start)<-colnames(Y)
    
    YN <- matrix(final.fit$resid,ncol=m,nrow=n) + matrix(c(0,final.fit$fit$fixef[-1])+final.fit$fit$fixef[1],nrow=n,ncol=m,byrow=T)
    dimnames(YN)<-dimnames(Y)
  }else{
    Sigma.cur=NULL
    YN <- matrix(fit2$resid,ncol=m,nrow=n) + matrix(c(0,fit2$fixef[-1])+fit2$fixef[1],nrow=n,ncol=m,byrow=T)
    dimnames(YN)<-dimnames(Y)
  }
  
  if(return.cov){
    list(Y=t(YN),Sigma=Sigma.cur)
  } else{
    t(YN)
  }
}


heatmap <- function(cor.mat){
  # Plots a red and blue heatmap
  # Args:
  #   cor.mat: a square matrix (e.g. a correlation matrix)
  # Returns:
  #   A plot is produced. No return arguments.
  
  f <- function(m) t(m)[,nrow(m):1]

  m <- nrow(cor.mat)
  
  ColorLevels <- seq(-1,1,length.out=101)
  ColorRamp <- rev(redblue(100))

  layout(matrix(c(1,1,2,0),2,2),widths=c(4,0.65))

  par(mar=c(3,1,1,3))
  image(1:m,1:m,f(cor.mat),col=ColorRamp,breaks=ColorLevels,xaxt = "n", yaxt = "n",ylab="",xlab="")
  axis(1, 1:m, labels = colnames(cor.mat), las = 2, line = -0.5, tick = 0)
  axis(4, 1:m, labels = rev(colnames(cor.mat)), las = 2, line = -0.5, tick = 0)

  par(mar = c(1,3,1.5,1))
  image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",cex.axis=1)
}
