library(igraph)
library(RSpectra)
library(mclust)
library(Matrix)
library(mvtnorm)

#source("emase.R")

reduceMix <- function(model)
{
   d <- model$d
   p <- model$parameters
   if(d==p$variance$d) return(model)
   p$mean <- model$parameters$mean[1:d,,drop=FALSE]
   p$variance$d <- d
   if(d==1){
      model$modelName <- 'V'
      p$variance$d <- 1
      p$variance$sigmasq <- p$variance$sigma[1,1,]
      p$variance$scale <- p$variance$sigmasq
      p$variance$modelName <- "V"
   } else {
      V <- array(0,dim=c(d,d,model$G))
      C <- array(0,dim=c(d,d,model$G))
      for(jj in 1:model$G){
         V[,,jj] <- p$variance$sigma[1:d,1:d,jj]
         C[,,jj] <- chol(V[,,jj])
      }
      p$variance$sigma <- V
      p$variance$cholsigma <- C
   }
   model$data <- model$data[,1:d,drop=FALSE]
   model$parameters <- p
   model
}

componentDensities <- function(x,pis,means,vars)
{
    a <- sapply(1:length(pis),function(i) {
          pis[i]*mvtnorm::dmvnorm(x,means[i,],vars[,,i])
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

## version of Mclust to run a bunch of times and return the best
mclust <- function(ntries=1,
                   adjustSubset=FALSE, # change the size of the subsets
                   maxSubset=8192,     # used for initialization
                   ...)
{
   if(adjustSubset)
      sub <- mclust.options("subset")
   m <- try(Mclust(...))
   if(ntries>1){
      best <- ifelse(inherits(m,'try-error'),-Inf,m$bic)
      if(adjustSubset){
         s <- sub+sample(0:max(0,maxSubset-sub),1)
         mclust.options(subset=s)
      }
      for(i in 2:ntries){
         z <- try(Mclust(...))
         if(!inherits(z,'try-error') && z$bic>best){
            m <- z
            best <- m$bic
         }
      }
   }
   if(adjustSubset)
      mclust.options(subset=sub)
   m
}

ase <- function(g,verbose=FALSE,adjust.diag=FALSE,normalize=FALSE,scale.by.values=FALSE,
                vector='u',max.d=12)
{
   if(verbose) cat("Extracting adjacency matrix\n")
   A <- as_adjacency_matrix(g)
   if(adjust.diag){
      if(verbose) cat("Adjusting the diagonal of the adjacency matrix\n")
      n <- gorder(g)
      A <- A + Diagonal(x=degree(g)/(n-1))
   }
   if(normalize){
      if(verbose) cat("Normalizing the adjacency matrix\n")
      d <- degree(g)
      dinv <- sqrt(1/d)
      dinv[d==0] <- 0
      Dinv <- Diagonal(x=dinv)
      A <- Dinv %*% A %*% Dinv
   }
   if(verbose) cat("Computing the spectrum\n")
   if(!is.directed(g)){
      sv <- eigs(A,k=max.d)
      z <- sv$vectors
      if(scale.by.values) {
         ev <- sv$values
         ev[ev<=0] <- 1
         z <- scale(z,center=FALSE,scale=1/sqrt(ev))
      }
   } else {
      sv <- svds(A,k=max.d)
      if(tolower(vector)=='u'){
         z <- sv$u
         if(scale.by.values) z <- scale(z,center=FALSE,scale=1/sqrt(sv$d))
      } else if(tolower(vector)=='v'){
         z <- sv$v
         if(scale.by.values) z <- scale(z,center=FALSE,scale=1/sqrt(sv$d))
      } else {
         error("unknown vector specification")
      }
   }
   z
}

aseCluster <- function(g,z,scale.by.values=FALSE,
                vector="u",     # only for directed grapohs
                only.embed=FALSE, # run ase()
                adjust.diag=FALSE,normalize=FALSE,
                min.d=1,
                min.k=1,
                max.k=12,max.d=max.k,
                short.circuit=FALSE, # if TRUE run Mclust to determing G
                #em=FALSE, # use the em code instead of the Mclust approx.
                maxIts=1000,epsilon=10^(-5),
                one.sigma=FALSE,   # this hasn't really been tested
                save.clusterings=TRUE,
                adjustSubset=FALSE,
                maxSubset=8192,
                au,av,    ## only for 'uv' and directed graphs
                ntries=1, ## number of tries of the Mclust algorithm
                verbose=FALSE){
   if(is.directed(g) && tolower(vector)=='uv' && only.embed==FALSE){
       ## Special case for 'uv'. This is not the only way to do this
       ## but it is the way that is currently implemented.
       if(missing(au)){
          au <- aseCluster(g,scale.by.values=scale.by.values,vector='u',
                   adjust.diag=adjust.diag,
                   normalize=normalize,
                   min.d=min.d,
                   min.k=min.k,
                   max.k=max.k,max.d=max.d,
                   short.circuit=short.circuit,
                   maxIts=maxIts,epsilon=epsilon,
                   one.sigma=one.sigma,
                   save.clusterings=save.clusterings,
                   ntries=ntries,
                   verbose=verbose)
       }
       if(missing(av)){
          av <- aseCluster(g,scale.by.values=scale.by.values,vector='v',
                   adjust.diag=adjust.diag,
                   normalize=normalize,
                   min.d=min.d,
                   min.k=min.k,
                   max.k=max.k,max.d=max.d,
                   short.circuit=short.circuit,
                   maxIts=maxIts,epsilon=epsilon,
                   one.sigma=one.sigma,
                   save.clusterings=save.clusterings,
                   ntries=ntries,
                   verbose=verbose)
      }
      z <- cbind(au$embedding[,1:au$d],av$embedding[,1:av$d])
      gs <- c(au$G,av$G)
      gs <- c(min(gs),sum(gs))
      m <- mclust(data=z,G=gs,verbose=FALSE,ntries=ntries,
                  adjustSubset=adjustSubset,
                  maxSubset=maxSubset)
      return(list(u=au,v=av,uv=m))
   }
   if(missing(z)){
      z <- ase(g,verbose=verbose,
               adjust.diag=adjust.diag,
               normalize=normalize,
               scale.by.values=scale.by.values,
                vector=vector,max.d=max.d)
   }
   if(only.embed) return(embedding)
   clusterings <- vector('list',max.d)
   ii <- 1
   BICS <- matrix(NA,nrow=max.d,ncol=max.k)
   NPARS <- matrix(NA,nrow=max.d,ncol=max.k)
   logn <- log(nrow(z))
   for(d in min.d:max.d){
      if(verbose) cat("d =",d,"\n")
      for(ncomp in min.k:max.k){
         if(verbose) cat("\tG =",ncomp,"\n")
         retries <- 0
         if(short.circuit){
            nc <- ncomp:max.k
         } else {
            nc <- ncomp
         }
         mix <- try(mclust(data=z[,1:d,drop=FALSE],G=nc,
                       verbose=verbose,
                       adjustSubset=adjustSubset,
                       maxSubset=maxSubset,
                       modelNames=ifelse(d==1,"V","VVV"),ntries=ntries))
         if(is.null(mix) || inherits(mix,'try-error')){
            warning("Could not fit the model.\n")
            # just fit a gaussian
            mix <- try(Mclust(data=z[,1:d,drop=FALSE],G=1,verbose=verbose,
                          modelNames=ifelse(d==1,"V","VVV")))
            if(is.null(mix) || inherits(mix,'try-error')){
               warning("Couldn't even get a one component solution.\n")
               mix <- list(modelName=ifelse(d==1,"V","VVV"),fit.failed=TRUE,
                           d=d,G=ncomp,bic=-Inf,loglik=NA,retries=retries)
            }
            mix$retries <- retries
            clusterings[[ii]] <- mix
            ii <- ii+1
         } else {
            mix$retries <- retries
            ncomp <- mix$G
            nm <- nMclustParams(mix$modelName,d=d,G=mix$G)
            if(d<ncol(z)){
               ## This is a bit of a kludge, but is probably good enough.
               ## It really should be recoded as in the other approach to
               ## fit a 0 mean spherical with fixed proportions.
               if(verbose) cat("fitting the noninformative part. d =",d,"\n")
               dat <- scale(z[,(d+1):ncol(z),drop=FALSE],scale=FALSE)
               mix$parameters$mean <- rbind(mix$parameters$mean,
                                            matrix(0,
                                                   nrow=ncol(dat),
                                                   ncol=mix$G))
               if(d == 1){
                  V <- mix$parameters$variance$sigmasq
               } else {
                  V <- mix$parameters$variance$sigma
               }
               VAR <- array(0,dim=c(max.d,max.d,mix$G))
               mix$parameters$variance <- list(modelName="VVV",
                    d=ncol(z),
                    G=mix$G,
                    sigma=VAR,
                    cholsigma=VAR)
               if(one.sigma){
                  s <- var(c(dat))
                  s <- rep(s,ncol(dat))
               }
               for(jj in 1:mix$G){
                  Z <- matrix(0,nrow=ncol(z),ncol=ncol(z))
                  if(d==1){
                     Z[1:d,1:d] <- V[jj]
                  } else {
                     Z[1:d,1:d] <- V[,,jj]
                  }
                  if(!one.sigma){
                     s <- sum(matrix(mix$z[,jj],nrow=1)%*%
                                     dat^2)/(ncol(dat)*sum(mix$z[,jj]))
                     s <- rep(s,ncol(dat))
                  } 
                  diag(Z) <- c(diag(Z[1:d,1:d,drop=FALSE]),s)
                  mix$parameters$variance$sigma[,,jj] <- Z
                  mix$parameters$variance$cholsigma[,,jj] <- chol(Z)
               }
               mix$modelName <- "VVV"
               mix$z <- mix$z[1,,drop=FALSE]
               mix$data <- z[1,,drop=FALSE]
               if(mix$G==1){
                  ll <- sum(dmvnorm(z,
                                    mean=mix$parameters$mean,
                                    mix$parameters$variance$sigma[,,1],
                                    log=TRUE))
               } else {
                  ll <- logLike(componentDensities(z,
                                                mix$parameters$pro,
                                                t(mix$parameters$mean),
                                                mix$parameters$variance$sigma))
               }
               mix$loglik <- ll
               bic <- 2*ll - (nm+ifelse(one.sigma,1,ncomp))*logn
               mix$bic <- bic
            } 
            if(verbose) cat("bic =",mix$bic,"\n")
            clusterings[[ii]] <- mix
            ii <- ii+1
            BICS[d,ncomp] <- mix$bic
            NPARS[d,ncomp] <- (nm+ifelse(d<ncol(z),
                                         ifelse(one.sigma,1,ncomp),0))
            if(verbose) {
               par(mfrow=c(2,2))
               image(1:max.d,1:max.k,BICS,ylab="G",
                     xlab="d",main="BIC",col=gray((255:0)/255))
               matplot(BICS,main="Components",ylab="BIC",xlab="d")
               matplot(t(BICS),main="Dimensions",ylab="BIC",xlab="G")
               image(1:max.d,1:max.k,NPARS,xlab="d",ylab="G",
                     main="#Params",col=gray((255:0)/255))
               par(mfrow=c(1,1))
            }
            if(short.circuit) break
         }
      }
   }
   b <- which.max(unlist(lapply(clusterings,'[[','bic')))
   model <- clusterings[[b]]
   if(is.infinite(model$bic)) {
      ## none of the models fit. Bail and try to fit anything
      ## on the full data set, any G, and pray.
      ## this should never happen.
      model <- try(mclust(data=z,G=1:max.k,verbose=verbose,ntries=ntries,
                          adjustSubset=adjustSubset,
                          maxSubset=maxSubset))
      if(inherits(model,'try-error'))
         error("Could not fit the model")
      clusterings <- NULL
   }
   if(verbose) cat("\tbic =",model$bic,"d =",model$d,"\n\tG =",model$G,"\n")
   membership1 <- predict(model,z)$classification
   membership2 <- try(predict(reduceMix(model),
                          z[,1:model$d,drop=FALSE])$classification)
   if(inherits(membership2,'try-error')) print(reduceMix(model))
   if(!save.clusterings) {
      clusterings <- NULL
   }
   list(d=model$d,
        call=match.call(),
        args=formals(aseCluster),
        n=gorder(g),
        s=gsize(g),
        embedding=z,
        clusterings=clusterings,
        membership=list(full=membership1,reduced=membership2),
        min.G=min.k,
        max.G=max.k,
        min.d=min.d,
        max.d=max.d,
        best.model=model,
        reduced.model=reduceMix(model),
        G=model$G,
        ntries=ntries,
        nparams=NPARS,
        best.index=b,
        bics=BICS)
}
