library(mclust)

firstPeak <- function(x) {
  n <- length(x)
  idx <- 1
  if (n > 1) {
    for (i in 2 : n) {
      if (x[i] - x[idx] < -6) {
        break
      }
      if (x[i] - x[idx] > 0) {
        idx <- i
      }
    }
  }
  return (idx)
}

# BIC value for GMM with 0 mean and equal diagonal covariance matrices 
# X: data matrices, n by d, each row is a sample
# k: number of components
bicGmmZero <- function(X, k) {
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  dd <- dim(X)[2]
  tao <- rep(1 / k, k)
  sigma2 <- rep(sum(diag(X %*% t(X)) / (dd * n)), k)
  factor <- 1:k
  sigma2 <- sigma2 * factor
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] / (sigma2[j]^(dd/2)) * exp(-sum(X[i,]^2) / (2 * sigma2[j]))
    }
  }
  loglik <- - (dd * n / 2) * log(2 * pi)
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  #print(loglik)
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      numerator <- 0
      for (i in 1 : n) {
        numerator <- numerator + t[i,j] * sum(X[i,]^2)
      }
      sigma2[j] <- numerator / (dd * sum(t[,j]))
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        f[i,j] <- tao[j] / (sigma2[j]^(dd/2)) * exp(-sum(X[i,]^2) / (2 * sigma2[j]))
      }
    }
    preLoglik <- loglik
    loglik <- - (dd * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  bic <- 2 * loglik - (2 * k - 1) * log(n)
  return (bic)
}

# use Jordan's model but iterate with inclusion and exclusion greedy algorithm
modelSelect2 <- function(X, kMax, dMax, bicDiff){
  maxBIC <- rep(-100, dMax)
  selectedIndex <- rep(0, dMax)
  selectedIndex[1] = 1
  
  # first index initialization
  for (k in 1:dMax)  {
    if (selectedIndex[k] == 0) {
      tempIndex <- selectedIndex
      tempIndex[k] <- 1
      mcObjSelected <- Mclust(X[, tempIndex == 1], 1:kMax)
      mcObjUnselected <- Mclust(X[, tempIndex == 0], 1:kMax)
      maxBIC[k] <- max(mcObjSelected$BIC, na.rm=T) + max(mcObjUnselected$BIC, na.rm=T)
    }
  }
  idx <- which.max(maxBIC)
  currentBIC <- max(maxBIC)
  selectedIndex[idx] <- 1
  
  terminate <- 0
  iter <- 0
  while (terminate == 0 && iter < 30) {
    terminate <- 1
    
    # inclusion step
    maxBIC <- rep(-100, dMax)
    for (k in 1:dMax) {
      if (selectedIndex[k] == 0) {
        tempIndex <- selectedIndex
        tempIndex[k] <- 1
        mcObjSelected <- Mclust(X[, tempIndex == 1], 1:kMax)
        mcObjUnselected <- Mclust(X[, tempIndex == 0], 1:kMax)
        maxBIC[k] <- max(mcObjSelected$BIC, na.rm=T) + max(mcObjUnselected$BIC, na.rm=T)
      }
    }
    idx <- which.max(maxBIC)
    if ((max(maxBIC) - currentBIC) > bicDiff) {
      selectedIndex[k] <- 1
      currentBIC = max(maxBIC)
      terminate <- 0
    }
    
    # exclusion step
    maxBIC <- rep(-100, dMax)
    for (k in 1:dMax) {
      if (selectedIndex[k] == 1) {
        tempIndex <- selectedIndex
        tempIndex[k] <- 0
        mcObjSelected <- Mclust(X[, tempIndex == 1], 1:kMax)
        mcObjUnselected <- Mclust(X[, tempIndex == 0], 1:kMax)
        maxBIC[k] <- max(mcObjSelected$BIC, na.rm=T) + max(mcObjUnselected$BIC, na.rm=T)
      }
    }
    idx <- which.max(maxBIC)
    if ((max(maxBIC) - currentBIC) > bicDiff) {
      selectedIndex[k] <- 0
      currentBIC = max(maxBIC)
      terminate <- 0
    }
    
    iter <- iter + 1
  }
  
  return(selectedIndex)
}

megaModelSelect <- function(X, G, dMax){
  #Do Mclust on selected variables only for X[,1:d], d = 1, 2, ... dMax
  mcObjs<-NULL
  
  #Do modeling of extra (not selected) variables. 
  idx <-1
  mcObjs[[1]] <- getMCforDim(X,idx,G)
  mcObjs[[1]] $BIC <- mcObjs[[1]]$BIC + getExtraBIC(X,idx)
  currPeak <- max(mcObjs[[1]]$BIC,na.rm=T)
  
  for(k in 2:dMax){
    mcObjs[[k]] <- getMCforDim(X,k,G)
    mcObjs[[k]]$BIC <- mcObjs[[k]]$BIC + getExtraBIC(X,k)
    propPeak = max(mcObjs[[k]]$BIC,na.rm=T)
    diff = propPeak - currPeak
    if(diff< -6){
      break
    } 
    if(diff > 0){
      idx <- k
      currPeak <- propPeak
    }
  }
  
  
  maxBICs <- lapply(mcObjs, function(mcObj.adj){
    max(mcObj.adj$BIC,na.rm=T)
  })
  
  dHat=idx
  return(list(dHat=dHat,
              kHat=mcObjs[[dHat]]$G,
              mcObjs=mcObjs, 
              maxBICs=maxBICs))
}

chooseDhat<-function(maxBICs){
  currPeak <- maxBICs[1]
  idx <-1
  for(k in 2:length(maxBICs)){
    diff = maxBICs[k] - currPeak
    if(diff< -5){
      return(idx)
    } 
    if(diff > 0){
      idx <- k
      currPeak <- maxBICs[k]
    }
  }
  return(idx)
}

getExtraBIC <- function(X,var){
  if(var==ncol(X)){
    mc.obj.fake<-NULL
    mc.obj.fake$BIC=0
    return(0)
  } else {
    max(Mclust(X[,(var+1):ncol(X)],G=1:10)$BIC,na.rm=T)
  }
}



getMCforDim <- function(X,d,G){
  mc.obj <- Mclust(X[,1:d], 1:G)
  return(mc.obj)
}

getBICGMM<-function(Xcol){
  return(max(Mclust(Xcol)$BIC,na.rm = T))
}
