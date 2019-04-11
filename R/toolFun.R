getElbows <- function(dat, n = 3, threshold = FALSE, plot = F) {
  if (is.matrix(dat))
    d <- sort(apply(dat,2,sd), decreasing=T)
  else
    d <- sort(dat,decreasing=T)
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) / (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < p) {
    q <- c(q, q + getElbows(d[(q+1):p], n=n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev")
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b")
      points(q,dat[q],col=2,pch=19)
    }
  }
  
  return(q)
}

nonpsd.laplacian <- function(A)
{
  n <- nrow(A)
  s <- rowSums(A)/(n-1)
  diag(A) <- diag(A)+s
  return(A)
}

stfp <- function(A, dim, scaling = TRUE)
{
  suppressMessages(require(irlba))
  
  if (dim(A)[1]==dim(A)[2]) { ## diagonal augmentation
    L <- nonpsd.laplacian(A)
  } else {
    L <- A
  }
  
  if(is.null(dim)) {
    L.svd <- irlba(L)
    #        L.svd <- svd(L)
    #        dim <- scree.thresh(L.svd$d) # ignore this, use "getElbows" instead
  } else {
    L.svd <- irlba(L,dim,dim)
    #        L.svd <- svd(L)
  }
  
  L.svd.values <- L.svd$d[1:dim]
  L.svd.vectors <- L.svd$v[,1:dim]
  
  if(scaling == TRUE) { # projecting to sphere
    if(dim == 1)
      L.coords <- sqrt(L.svd.values) * L.svd.vectors
    else
      L.coords <- L.svd.vectors %*% diag(sqrt(L.svd.values))
  }
  else {
    L.coords <- L.svd.vectors
  }
  
  return(list(X=L.coords,D=L.svd.values))
}

eigASE <- function(A, dim, scaling = TRUE)
{
  if (dim(A)[1]==dim(A)[2]) { ## diagonal augmentation
    L <- nonpsd.laplacian(A)
  } else {
    L <- A
  }
  L.eig <- eigen(L)
  L.eig.values <- L.eig$values[1:dim]
  L.eig.vectors <- L.eig$vectors[,1:dim]
  
  if(scaling == TRUE) {
    if(dim == 1)
      L.coords <- sqrt(L.eig.values) * L.eig.vectors
    else
      L.coords <- L.eig.vectors %*% diag(sqrt(L.eig.values))
  }
  else {
    L.coords <- L.eig.vectors
  }
  
  return(list(X=L.coords,D=L.eig.values))
}

oos.ase <- function(A11, A12, d)
{
  require(irlba)
  A11.svd <- irlba(A11,d)
  Xhat.insample <- A11.svd$u %*% diag(sqrt(A11.svd$d))
  Xhat.oos <- A12 %*% A11.svd$u%*%diag(1/sqrt(A11.svd$d))
  return(list(Xhat.insample=Xhat.insample, Xhat.oos=Xhat.oos))
}