library(igraph)
library(RSpectra)
library(mclust)
library(Matrix)
library(mvtnorm)

reduceMixEM <- function(mix)
{
   d <- mix$d
   return(list(pis=mix$pis,means=mix$means[,1:d],vars=mix$vars[1:d,1:d,],
               d=d))
}


componentDensities <- function(x,pis,means,vars,d)
{
    if(missing(d)) d <- ncol(x)
    a <- sapply(1:length(pis),function(i) {
          pis[i]*mvtnorm::dmvnorm(x[,1:d],means[i,1:d],vars[1:d,1:d,i])
       })
    a[is.nan(a)] <- 0
    a
}

logLike <- function(cd,min.likelihood=10^(-10))
{
      a <- apply(cd,1,sum)
      a[a<min.likelihood] <- min.likelihood
      sum(log(a))
}

clusterMix <- function(x,mix,d) 
{
   if(missing(d)) d <- ncol(x)
   apply(componentDensities(x[,1:d,drop=FALSE],mix$pis,mix$means,mix$vars,d=d),
         1,which.max)
}

taus <- function(cd)
{
   apply(cd,1,function(x) x/sum(x))
}

means <- function(x,tau)
{
   t(sapply(1:nrow(tau),function(i) apply(tau[i,]*x,2,sum)/sum(tau[i,])))
}

## #params = #pis - 1 + k*(mean + informative cov + non-informative)
## in the case of one sigma for all groups, there is only one parameter for the non-informative
nParams <- function(ncomps,d,max.d,one.sigma=FALSE){
   nf <- as.integer(d<max.d && !one.sigma)
   #p           G    (mean   var       non-inf sigma^2)
   ncomps-1 + ncomps*(d + d*(d+1)/2+nf) + as.integer(one.sigma)
}

vars <- function(x,m,d,tau,one.sigma)
{
   n <- nrow(x)
   dd <- ncol(x)
   v <- array(0,dim=c(dd,dd,nrow(tau)))
   for(i in 1:nrow(tau)){
      y <- scale(x,center=m[i,],scale=FALSE)
      if(d<dd){
         dat <- scale(x[,(d+1):dd,drop=FALSE],scale=FALSE)
      }
      for(j in 1:n){
         v[1:d,1:d,i] <- v[1:d,1:d,i] + tau[i,j]*t(y[j,1:d,drop=FALSE])%*%
                                                 y[j,1:d,drop=FALSE]
      }
      v[,,i] <- v[,,i]/sum(tau[i,])
      if(d<dd){
         if(one.sigma){
           s <- var(c(dat))
         } else {
            s <- sum(tau[i,]%*%dat^2)/(ncol(dat)*sum(tau[i,]))
            s <- rep(s,ncol(dat))
         }
         if(d==1) s <- c(v[1,1,i],s)
         else s <- c(diag(v[1:d,1:d,i]),s)
         diag(v[,,i]) <- s
      }
   }
   v
}

fitMix <- function(data,d,ncomps,one.sigma=FALSE,maxIts=1000,epsilon=10^(-5),verbose=FALSE,tries=3)
{
   ## Use Mclust to get a starting point for the full mixture model.
   if(verbose) cat("Initializing the informative parameters\n")
   m <- try(Mclust(data[,1:d,drop=FALSE],G=ncomps,modelName=ifelse(d==1,'V','VVV'),verbose=verbose),silent=TRUE)
   ## Sometimes Mclust fails to find a model. Retry it a couple times until it does, or bail
   ## out and do something naive to initialize.
   if(inherits(m,'try-error') || is.null(m)){
      if(verbose) cat("trying again\n")
      G <- ncomps
      for(j in 1:tries){
         inds <- sample(nrow(data),nrow(data),replace=TRUE)
         x <- data[inds,1:d,drop=FALSE]
         if(j>1){
            ## dither the points a bit
            x <- x+rnorm(nrow(data)*d,0,0.001)
         }
         ## try mclust again
         m <- try(Mclust(x,G=ncomps,modelName=ifelse(d==1,'V','VVV'),verbose=verbose),silent=TRUE)
         if(!inherits(m,'try-error')) break
      }
      ## if it still fails, bail and do something arbitrary
      if(inherits(m,'try-error') || is.null(m)){
         if(verbose) cat("Bailing -- could not initialize with Mclust\n")
         ps <- rep(1/G,G)
         ms <- data[sample(nrow(data),G),1:d]
         vs <- array(0,dim=c(d,d,G))
         s <- sd(c(data))/G
         for(j in 1:G)
            diag(vs[,,j]) <- s
      } else {
         ms <- t(m$parameters$mean)
         vs <- m$parameters$variance$sigma
         ps <- m$parameters$pro
      }
   } else {
      ms <- t(m$parameters$mean)
      vs <- m$parameters$variance$sigma
      ps <- m$parameters$pro
   }
   ## add 0 mean, diag covariance to the model for the uninformative dimensions
   if(d<ncol(data)){
      ## initialize the uninformative variances to the overall variance
      ## set the uninformative means to 0
      if(verbose) cat("Initializing the non-informative parameters\n")
      tmp <- matrix(0,ncol=ncol(data),nrow=ncomps)
      tmp[,1:d] <- ms
      ms <- tmp
      tmp <- array(0,dim=c(ncol(data),ncol(data),ncomps))
      v <- var(c(data[,-(1:d)]))
      for(i in 1:ncomps){
         if(d==1) {
            tmp[1,1,i] <- vs[i]
         } else {
            tmp[1:d,1:d,i] <- vs[,,i]
         }
         if(d<(ncol(data)-1)){
            diag(tmp[(d+1):ncol(data),(d+1):ncol(data),i]) <- v
         } else {
            tmp[ncol(data),ncol(data),i] <- v
         }
      }
      vs <- tmp
   } 
   logn <- log(nrow(data))
   cd <- componentDensities(data,ps,ms,vs)
   ll <- logLike(cd)
   if(is.na(ll)) error("invalid initial log likelihood")
   for(i in 1:maxIts){
      tau <- taus(cd)
      nm <- matrix(0,nrow=nrow(ms),ncol=ncol(ms))
      nm[,1:d] <- means(data[,1:d,drop=FALSE],tau)
      nv <- vars(data,ms,d,tau,one.sigma)
      ps <- apply(tau,1,sum)/nrow(data)
      cd <- componentDensities(data,ps,nm,nv)
      ll2 <- logLike(cd)
      ms <- nm
      vs <- nv
      if(verbose) cat("\t",i,ll2,ll,ll2-ll,"\n")
      if(is.na(ll2)) error("invalid log likelihood")
      if((ll2-ll)<epsilon) break
      ll <- ll2
   }
   nf <- as.integer(d<ncol(data))
   bic <- 2*ll2 - nParams(ncomps,d,ncol(data),one.sigma=one.sigma)*logn
   list(pis=ps,means=ms,vars=vs,d=d,G=ncomps,bic=bic,likelihood=ll)
}

aseClusterEM <- function(g,z,scale.by.values=TRUE,vectors="u",
                only.embed=FALSE,
                adjust.diag=FALSE,normalize=FALSE,
                remove.zeros=FALSE,
                zero.tol=10^-6,
                min.d=1,
                min.k=1,
                max.k=12,max.d=max.k,
                one.sigma=FALSE,
                maxIts=1000,epsilon=10^(-5),
                max.tries=5,
                save.embedding=TRUE,
                bail=TRUE, # bail out if the bic decreases
                verbose=FALSE){
   if(missing(z)){
      z <- ase(g,verbose=verbose,
               adjust.diag=adjust.diag,
               normalize=normalize,
               scale.by.values=scale.by.values,
                vector=vector,max.d=max.d)
   }
   if(missing(z)){
      z <- ase(g,verbose=verbose,
               adjust.diag=adjust.diag,
               normalize=normalize,
               scale.by.values=scale.by.values,
                vector=vector,max.d=max.d)
   }
   if(remove.zeros){
      zinds <- which(apply(z,1,function(x) sum(x^2)<zero.tol))
      if(length(zinds)>0) z <- z[-zinds,]
   } else {
      zinds <- NULL
   }
   if(only.embed) return(list(embedding=z,zero.indices=zinds))
   clusterings <- vector('list',max.d)
   ii <- 1
   BICS <- matrix(NA,nrow=max.d,ncol=max.k)
   NPARS <- matrix(NA,nrow=max.d,ncol=max.k)
   logn <- log(nrow(z))
   for(d in min.d:max.d){
      if(verbose) cat("d =",d,"\n")
      if(min.k==1){
         ## First do the G=1 case of a single Gaussian
         if(verbose) cat("G = 1\n")
         m <- apply(z,2,mean)
         v <- var(z)
         if(d<max.d) {
            m[-(1:d)] <- 0
            v2 <- var(c(z[,(d+1):max.d]))
            v[1:d,-(1:d)] <- 0
            v[-(1:d),1:d] <- 0
            v[-(1:d),-(1:d)] <- 0
            if((d+1)==max.d){
               v[max.d,max.d] <- v2
            } else {
               diag(v[(d+1):max.d,(d+1):max.d]) <- v2
            }
         }
         ll2 <- sum(log(mvtnorm::dmvnorm(z,m,v)))
         nf <- as.integer(d<ncol(z))
         #bic <- 2*ll2-(d+d*(d+1)/2+nf)*logn
         bic <- 2*ll2 - nParams(1,d,max.d,one.sigma=one.sigma)*logn
         if(verbose) cat("bic =",bic,"\n")
         clusterings[[ii]] <- list(pis=1,means=m,vars=v,d=d,G=1,bic=bic)
         ii <- ii+1
         BICS[d,1] <- bic
         NPARS[d,1] <- nParams(1,d,max.d,one.sigma=one.sigma)
         min.k <- 2
      }
      for(ncomps in min.k:max.k){
         if(verbose) cat("G =",ncomps,"\n")
         mix <- try(fitMix(z,d,ncomps,one.sigma=one.sigma,maxIts=maxIts,epsilon=epsilon,verbose=verbose),
                    silent=TRUE)
         iii <- 1
         while(inherits(mix,'try-error')){
            ## Try again!
            inds <- sample(nrow(z),nrow(z),replace=TRUE)
            x <- z[inds,]
            ## dither the points
            x <- x+rnorm(prod(dim(z)),0,0.001)
            mix <- try(fitMix(x,d,ncomps,
                              one.sigma=one.sigma,maxIts=maxIts,epsilon=epsilon,verbose=verbose),silent=TRUE)
            iii <- iii+1
            if(iii>max.tries) break
         }
         if(!inherits(mix,'try-error')){
            if(verbose) cat("bic =",mix$bic,"\n")
            clusterings[[ii]] <- mix
            ii <- ii+1
            BICS[d,ncomps] <- mix$bic
            NPARS[d,ncomps] <- nParams(ncomps,d,max.d,one.sigma=one.sigma)
            if(verbose) {
               par(mfrow=c(2,2))
               image(1:max.d,1:max.k,BICS,ylab="G",xlab="d",main="BIC",col=gray((255:0)/255))
               matplot(BICS,main="Components",ylab="BIC",xlab="d")
               matplot(t(BICS),main="Dimensions",ylab="BIC",xlab="G")
               image(1:max.d,1:max.k,NPARS,xlab="d",ylab="G",main="#Params",col=gray((255:0)/255))
               par(mfrow=c(1,1))
            }
            if(bail && !is.na(BICS[d,ncomps-1]) && (mix$bic<BICS[d,ncomps-1])){
               break
            }
         }
      }
   }
   b <- which.max(unlist(lapply(clusterings,'[[','bic')))
   model <- clusterings[[b]]
   if(verbose) cat("\tbic =",model$bic,"d =",model$d,"\n\tG =",model$G,"\n")
   membership1 <- clusterMix(z,model)
   membership2 <- clusterMix(z[,1:model$d],reduceMixEM(model))
   if(!save.embedding) z <- NULL
   list(d=model$d,
        call=match.call(),
        args=formals(aseCluster),
        n=gorder(g),
        s=gsize(g),
        embedding=z,
        zero.indices=zinds,
        clusterings=clusterings,
        membership=list(full=membership1,reduced=membership2),
        min.G=min.k,
        max.G=max.k,
        min.d=min.d,
        max.d=max.d,
        best.model=model,
        reduced.model=reduceMixEM(model),
        G=model$G,
        nparams=NPARS,
        best.index=b,
        bics=BICS)
}
