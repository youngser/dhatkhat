# Selecting (d, K) pair by the maximum BIC
# can choose global or local maximum, globalMax = 0 means it is the global maximum, otherwise it means
# the index of the local maximum from the left in terms of d
# bic is a dd * kmax matrix, whose d,k-th entry is the BIC value for each (d, k) pair
maxBIC1 <- function(bic, globalMax = 0, regression = 0, penalty = 1, n = 0, dmax = 0) {
  dd <- dim(bic)[1]
  kmax <- dim(bic)[2]
  dhat <- 1
  khat <- 1
  
  if (penalty > 1) {
    for (d in 1 : dd) {
      for (k in 1 : kmax) {
        numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2 + dmax - d)
        bic[d,k] <- bic[d,k] - (penalty - 1) * numOfParameters * log(n)
      }
    }
  }
  
  if ((globalMax == 0) && (regression == 0)) {
    maxBic <- 1e-10
    for (d in 1 : dd) {
      for (k in 1 : kmax) {
        if (bic[d, k] > maxBic) {
          maxBic <- bic[d, k]
          dhat <- d
          khat <- k
        }
      }
    }
    return (list(dhat = dhat,
                 khat = khat,
                 maxBic = maxBic))
  } else if ((globalMax == 0) && (regression > 0)) {
    khats <- rep(0, dd)
    maxBics <- rep(1e-10, dd)
    for (d in 1 : dd) {
      for (k in 1 : kmax) {
        if (bic[d, k] > maxBics[d]) {
          maxBics[d] <- bic[d, k]
          khats[d] <- k
        }
      }
    }
    x <- 1:dd
    y <- t(t(maxBics))
    # X <- matrix(rep(0, dd*3), dd, 3)
    # X[,1] <- x^2
    # X[,2] <- x
    # X[,3] <- rep(1, dd)
    X <- matrix(rep(0, dd* (regression + 1)), dd, (regression + 1))
    for (i in 1 : (regression + 1)) {
      X[,i] <- x^(regression - i + 1)
    }
    beta <- solve(t(X) %*% X) %*% t(X) %*% y
    # dhat <- round(- beta[2] / (2 * beta[1]))
    # if (beta[1] < 0) {
    #   if (dhat < 1) {
    #     dhat = 1
    #   } else if (dhat > dd) {
    #     dhat = dd
    #   }
    # } else if (beta[2] > 0) {
    #   if (maxBics[dd] > maxBics[1]) {
    #     dhat = dd
    #   } else {
    #     dhat = 1
    #   }
    # }
    dhat <- which.max(X %*% beta)
    khat <- khats[dhat]
    return (list(dhat = dhat,
                 khat = khat,
                 maxBics = maxBics,
                 beta = beta))
  } else if (globalMax > 0) {
    khats <- rep(0, dd)
    isPeak <- rep(0, dd)
    maxBics <- rep(1e-10, dd)
    for (d in 1 : dd) {
      for (k in 1 : kmax) {
        if (bic[d, k] > maxBics[d]) {
          maxBics[d] <- bic[d, k]
          khats[d] <- k
        }
      }
    }
    if (dd > 2) {
      for (i in 2 : (dd - 1)) {
        if ((maxBics[i] > maxBics[i-1]) && (maxBics[i] > maxBics[i+1])) {
          isPeak[i] <- 1
        }
      }
      if (maxBics[dd] > maxBics[dd-1]) {
        isPeak[dd] <- 1
      }
    } else if (dd == 2) {
      if (maxBics[2] > maxBics[1]) {
        isPeak[2] <- 1
      } else {
        isPeak[1] <- 1
      }
    } else {
      isPeak[1] <- 1
    }
    
    index <- 1
    peaks <- rep(0, sum(isPeak))
    for (i in 1 : dd) {
      if (isPeak[i] == 1) {
        peaks[index] <- i
        index <- index + 1
      }
    }
    if (globalMax > length(peaks)) {
      globalMax <- length(peaks)
    }
    return (list(dhat = peaks[globalMax],
                 khat = khats[peaks[globalMax]],
                 maxBics = maxBics))
  }
}

# calulate the BIC of the blockwise model 1, and then do clustering
# \mu_{d+1} = ... = \mu_{dmax} = 0; \sigma_{d+1}, ... , \sigma_{dmax} don't need to be equal
bicGmmBlock1 <- function(X, d, k, tao = 0) {
  # print("Block1")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  dmax <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu <- NULL
  Sigma1 <- NULL
  sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), 1 : d, drop = FALSE]
    X2 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), (d + 1) : dmax, drop = FALSE]
    head <- head + tao[i]
    mu[[i]] <- t(t(colMeans(X1)))
    Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
    for (l in 1 : dim(X1)[1]) {
      Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*% 
        (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
    }
    sigma2[[i]] <- rep(0, dmax - d)
    for (j in 1 : (dmax - d)) {
      for (l in 1 : dim(X2)[1]) {
        sigma2[[i]][j] <- sigma2[[i]][j] + X2[l, j]^2 / dim(X2)[1]
      }
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * prod(sigma2[[j]])^(-0.5) * 
        exp(-(X[i, 1 : d, drop = FALSE] - t(mu[[j]])) %*% 
              solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) / 2 - 
              X[i, (d + 1) : dmax, drop = FALSE] %*% diag(1 / sigma2[[j]], dmax - d, dmax - d) %*%
              t(X[i, (d + 1) : dmax, drop = FALSE]) / 2)
    }
  }
  loglik <- -(dmax * n / 2) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  #print(loglik)
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu[[j]] <- matrix(rep(0, d), d, 1)
      for (i in 1 : n) {
        mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, 1 : d, drop = FALSE])
      }
      mu[[j]] <- mu[[j]] / sum(t[,j])
      Sigma1[[j]] <- matrix(rep(0, d^2), d, d)
      sigma2[[j]] <- rep(0, dmax - d)
      for (i in 1 : n) {
        Sigma1[[j]] <- Sigma1[[j]] + t[i,j] * (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) %*% 
          (X[i, 1 : d, drop = FALSE] - t(mu[[j]]))
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * X[i, (d + 1) : dmax]^2
      }
      Sigma1[[j]] <- Sigma1[[j]] / sum(t[,j])
      sigma2[[j]] <- sigma2[[j]] / sum(t[,j])
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * prod(sigma2[[j]])^(-0.5) * 
          exp(-(X[i, 1 : d, drop = FALSE] - t(mu[[j]])) %*% 
                solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) / 2 - 
                X[i, (d + 1) : dmax, drop = FALSE] %*% diag(1 / sigma2[[j]], dmax - d, dmax - d) %*%
                t(X[i, (d + 1) : dmax, drop = FALSE]) / 2)
      }
    }
    preLoglik <- loglik
    loglik <- - (dmax * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2 + dmax - d)
  bic <- 2 * loglik - numOfParameters * log(n)
  classification <- rep(0, n)
  for (i in 1 : n) {
    classification[i] <- which.max(f[i,])
  }
  return (list(bic = bic,
               classification = classification,
               d = d,
               k = k))
}

# calulate the BIC of the blockwise model 2, and then do clustering
# \mu_{d+1} = ... = \mu_{dmax} = 0; \sigma_{d+1} = ... = \sigma_{dmax}
# at least of init and k needs to be assigned
bicGmmBlock2 <- function(X, d, k, init = NULL) {
  # print("Block2")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  dmax <- dim(X)[2]
  
  # initialization of parameters
  if (is.null(init)) {
    tao <- rep(1 / k, k)
    mu <- NULL
    Sigma1 <- NULL
    sigma2 <- NULL
    head <- 0
    for (i in 1 : k) { 
      X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), 1 : d, drop = FALSE]
      X2 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), (d + 1) : dmax, drop = FALSE]
      head <- head + tao[i]
      mu[[i]] <- t(t(colMeans(X1)))
      Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
      for (l in 1 : dim(X1)[1]) {
        Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*% 
          (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
      }
      sigma2[[i]] <- 0
      for (l in 1 : dim(X2)[1]) {
        sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
      }
    }
  } else {
    label <- unique(init)
    numOfPoints <- rep(0, k)
    for (i in 1 : length(label)) {
      numOfPoints[i] <- sum(init == label[i])
    }
    indices <- sort.list(numOfPoints, decreasing = TRUE)
    if (length(label) >= k) {
      tao <- rep(0, k)
      mu <- NULL
      Sigma1 <- NULL
      sigma2 <- NULL
      for (i in 1 : k) {
        tao[i] <- sum(init == label[indices[i]])
        X1 <- X[init == label[indices[i]], 1 : d, drop = FALSE]
        X2 <- X[init == label[indices[i]], (d + 1) : dmax, drop = FALSE]
        mu[[i]] <- t(t(colMeans(X1)))
        Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
        for (l in 1 : dim(X1)[1]) {
          Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*% 
            (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
        }
        sigma2[[i]] <- 0
        for (l in 1 : dim(X2)[1]) {
          sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
        }
      }
      tao <- tao / sum(tao)
    } else {
      # length(label) < k
      tao <- rep(0, k)
      mu <- NULL
      Sigma1 <- NULL
      sigma2 <- NULL
      XX <- X[init == label[indices[1]], , drop = FALSE]
      numOfParti <- k - length(label) + 1
      numOfRows <- sum(init == label[indices[1]])
      for (i in 1 : numOfParti) {
        head <- floor(numOfRows * (i - 1) / numOfParti) + 1
        tail <- floor(numOfRows * (i) / numOfParti)
        tao[i] <- tail - head + 1
        X1 <- XX[head : tail, 1 : d, drop = FALSE]
        X2 <- XX[head : tail, (d + 1) : dmax, drop = FALSE]
        mu[[i]] <- t(t(colMeans(X1)))
        Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
        for (l in 1 : dim(X1)[1]) {
          Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*% 
            (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
        }
        sigma2[[i]] <- 0
        for (l in 1 : dim(X2)[1]) {
          sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
        }
      }
      for (i in (numOfParti + 1) : k) {
        tao[i] <- sum(init == label[indices[i - numOfParti + 1]])
        X1 <- X[init == label[indices[i - numOfParti + 1]], 1 : d, drop = FALSE]
        X2 <- X[init == label[indices[i - numOfParti + 1]], (d + 1) : dmax, drop = FALSE]
        mu[[i]] <- t(t(colMeans(X1)))
        Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
        for (l in 1 : dim(X1)[1]) {
          Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*% 
            (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
        }
        sigma2[[i]] <- 0
        for (l in 1 : dim(X2)[1]) {
          sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
        }
      }
      tao <- tao / sum(tao)
    }
  }
  
  # Check if there is singularity
  for (j in 1 : k) {
    if (kappa(Sigma1[[j]]) > 1e10) {
      print("Singular at initialization")
      return (list(bic = 0,
                   classification = init,
                   d = d,
                   k = k,
                   iter = 0,
                   logliks = NULL,
                   tao = tao,
                   mu = mu,
                   Sigma1 = Sigma1,
                   sigma2 = sigma2))
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  fTemp <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * sigma2[[j]]^(-(dmax - d) / 2) * 
        exp(-(X[i, 1 : d, drop = FALSE] - t(mu[[j]])) %*% 
              solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) / 2 - 
              sum(X[i, (d + 1) : dmax]^2) / (2 * sigma2[[j]]))
    }
  }
  loglik <- -(dmax * n / 2) * log(2 * pi)
  for (i in 1 : n) {
    fi <- sum(f[i,])
    if (fi > 1e-1000) {
      loglik <- loglik + log(fi)
    } else {
      loglik <- loglik - 1000
    }
  }
  #print(loglik)
  
  flag <- 0
  logliks <- rep(0, 100)
  logliks[iter + 1] <- loglik
  while ((flag == 0) && (diff > 1e-5) && (iter < 100)) {
    # print(iter)
    # if (iter == 8) {
    #   browser()
    # }
    # E step
    for (i in 1 : n) {
      if (sum(f[i,]) == 0) {
        for (j in 1 : k) {
          t[i,j] <- 1/k
        }
      } else {
        for (j in 1 : k) {
          t[i,j] <- f[i,j] / sum(f[i,])
          if (t[i,j] < 1e-10) {
            t[i,j] <- 1e-10
          }
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu[[j]] <- matrix(rep(0, d), d, 1)
      for (i in 1 : n) {
        mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, 1 : d, drop = FALSE])
      }
      mu[[j]] <- mu[[j]] / sum(t[,j])
      Sigma1[[j]] <- matrix(rep(0, d^2), d, d)
      sigma2[[j]] <- 0
      for (i in 1 : n) {
        Sigma1[[j]] <- Sigma1[[j]] + t[i,j] * (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) %*% 
          (X[i, 1 : d, drop = FALSE] - t(mu[[j]]))
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * sum(X[i, (d + 1) : dmax]^2)
      }
      Sigma1[[j]] <- Sigma1[[j]] / sum(t[,j])
      sigma2[[j]] <- sigma2[[j]] / ((dmax - d) * sum(t[,j]))
      # Calculate determinant and inverse of Sigma1 (181104)
    }
    # check convergence
    detSigma1 <- NULL
    invSigma1 <- NULL
    for (j in 1 : k) {
      if (kappa(Sigma1[[j]]) > 1e10) {
        print("Singular in the EM steps")
        numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2 + 1)
        #loglik <- preLoglik
        bic <- 2 * loglik - numOfParameters * log(n)
        classification <- rep(0, n)
        for (i in 1 : n) {
          classification[i] <- which.max(f[i,])
        }
        return (list(bic = bic,
                     classification = classification,
                     d = d,
                     k = k,
                     iter = iter,
                     logliks = logliks,
                     tao = tao,
                     mu = mu,
                     Sigma1 = Sigma1,
                     sigma2 = sigma2))
      }
      detSigma1[[j]] <- myDet(Sigma1[[j]])
      invSigma1[[j]] <- solve(Sigma1[[j]])
    }
    
    for (i in 1 : n) {
      for (j in 1 : k) {
        if ((abs(detSigma1[[j]]) < 1e-1000) || (abs(sigma2[[j]]) < 1e-1000)) {
          flag <- 1
        } else {
          fTemp[i,j] <- tao[j] * ((detSigma1[[j]])^(-0.5)) * sigma2[[j]]^(-(dmax - d) / 2) * 
            exp(-(X[i, 1 : d, drop = FALSE] - t(mu[[j]])) %*% 
                  invSigma1[[j]] %*% (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) / 2 - 
                  sum(X[i, (d + 1) : dmax]^2) / (2 * sigma2[[j]]))
        }
      }
    }
    if (flag == 0) {
      f <- fTemp
    }
    preLoglik <- loglik
    loglik <- - (dmax * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      fi <- sum(f[i,])
      if (fi > 1e-1000) {
        loglik <- loglik + log(fi)
      } else {
        loglik <- loglik - 1000
      }
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
    logliks[iter + 1] <- loglik
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2 + 1)
  bic <- 2 * loglik - numOfParameters * log(n)
  classification <- rep(0, n)
  for (i in 1 : n) {
    classification[i] <- which.max(f[i,])
  }
  return (list(bic = bic,
               classification = classification,
               d = d,
               k = k,
               iter = iter,
               logliks = logliks,
               tao = tao,
               mu = mu,
               Sigma1 = Sigma1,
               sigma2 = sigma2))
}

# compute the posterior density of a given vector on one component of GMM with GMMBlock2 pattern
densityBlock2 <- function(x, tao, mu, Sigma, sigma2) {
  dmax <- length(x)
  d <- dim(Sigma)[1]
  if (kappa(Sigma) > 1e10) {
    f <- 0
  } else {
    f <- tao * (2 * pi)^(-dmax/2) * (myDet(Sigma)^(-0.5)) * sigma2^(-(dmax - d) / 2) * 
      exp(-(t(x[1 : d]) - t(mu)) %*% 
            solve(Sigma) %*% (t(t(x[1 : d])) - t(t(mu))) / 2 - 
            sum(x[(d + 1) : dmax]^2) / (2 * sigma2))
  }
  return (f)
}

# calulate the BIC of the blockwise model 3, and then do clustering
# \mu_{d+1}, ... ,\mu_{dmax} have no constraints; \sigma_{d+1}, ... , \sigma_{dmax} don't need to
# be equal
bicGmmBlock3 <- function(X, d, k, tao = 0) {
  # print("Block3")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  dmax <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu1 <- NULL
  mu2 <- NULL
  Sigma1 <- NULL
  sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), 1 : d, drop = FALSE]
    X2 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), (d + 1) : dmax, drop = FALSE]
    head <- head + tao[i]
    mu1[[i]] <- t(t(colMeans(X1)))
    mu2[[i]] <- t(t(colMeans(X2)))
    Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
    for (l in 1 : dim(X1)[1]) {
      Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu1[[i]]) %*% 
        (X1[l, , drop = FALSE] - t(mu1[[i]])) / dim(X1)[1]
    }
    sigma2[[i]] <- rep(0, dmax - d)
    for (j in 1 : (dmax - d)) {
      for (l in 1 : dim(X2)[1]) {
        sigma2[[i]][j] <- sigma2[[i]][j] + (X2[l, j] - mu2[[i]][j])^2 / dim(X2)[1]
      }
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * prod(sigma2[[j]])^(-0.5) * 
        exp(-(X[i, 1 : d, drop = FALSE] - t(mu1[[j]])) %*% 
              solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu1[[j]]) / 2 - 
              (X[i, (d + 1) : dmax, drop = FALSE] - t(mu2[[j]])) %*%
              diag(1 / sigma2[[j]], dmax - d, dmax - d) %*%
              (t(X[i, (d + 1) : dmax, drop = FALSE]) - mu2[[j]]) / 2)
    }
  }
  loglik <- -(dmax * n / 2) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  #print(loglik)
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu1[[j]] <- matrix(rep(0, d), d, 1)
      mu2[[j]] <- matrix(rep(0, dmax - d), (dmax - d), 1)
      for (i in 1 : n) {
        mu1[[j]] <- mu1[[j]] + t[i,j] * t(X[i, 1 : d, drop = FALSE])
        mu2[[j]] <- mu2[[j]] + t[i,j] * t(X[i, (d + 1) : dmax, drop = FALSE])
      }
      mu1[[j]] <- mu1[[j]] / sum(t[,j])
      mu2[[j]] <- mu2[[j]] / sum(t[,j])
      Sigma1[[j]] <- matrix(rep(0, d^2), d, d)
      sigma2[[j]] <- rep(0, dmax - d)
      for (i in 1 : n) {
        Sigma1[[j]] <- Sigma1[[j]] + t[i,j] * (t(X[i, 1 : d, drop = FALSE]) - mu1[[j]]) %*% 
          (X[i, 1 : d, drop = FALSE] - t(mu1[[j]]))
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * (X[i, (d + 1) : dmax] - t(mu2[[j]])[1, ])^2
      }
      Sigma1[[j]] <- Sigma1[[j]] / sum(t[,j])
      sigma2[[j]] <- sigma2[[j]] / sum(t[,j])
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * prod(sigma2[[j]])^(-0.5) * 
          exp(-(X[i, 1 : d, drop = FALSE] - t(mu1[[j]])) %*% 
                solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu1[[j]]) / 2 - 
                (X[i, (d + 1) : dmax, drop = FALSE] - t(mu2[[j]])) %*%
                diag(1 / sigma2[[j]], dmax - d, dmax - d) %*%
                (t(X[i, (d + 1) : dmax, drop = FALSE]) - mu2[[j]]) / 2)
      }
    }
    preLoglik <- loglik
    loglik <- - (dmax * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (dmax + (d + 1) * d / 2 + dmax - d)
  bic <- 2 * loglik - numOfParameters * log(n)
  classification <- rep(0, n)
  for (i in 1 : n) {
    classification[i] <- which.max(f[i,])
  }
  return (list(bic = bic,
               classification = classification,
               d = d,
               k = k))
}

# calulate the BIC of the blockwise model 4, and then do clustering
# \mu_{d+1}, ... ,\mu_{dmax} have no constraints; \sigma_{d+1} = ... = \sigma_{dmax}
bicGmmBlock4 <- function(X, d, k, tao = 0) {
  # print("Block4")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  dmax <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu1 <- NULL
  mu2 <- NULL
  Sigma1 <- NULL
  sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), 1 : d, drop = FALSE]
    X2 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), (d + 1) : dmax, drop = FALSE]
    head <- head + tao[i]
    mu1[[i]] <- t(t(colMeans(X1)))
    mu2[[i]] <- t(t(colMeans(X2)))
    Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
    for (l in 1 : dim(X1)[1]) {
      Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu1[[i]]) %*% 
        (X1[l, , drop = FALSE] - t(mu1[[i]])) / dim(X1)[1]
    }
    sigma2[[i]] <- 0
    for (l in 1 : dim(X2)[1]) {
      sigma2[[i]] <- sigma2[[i]] + sum((X2[l, ] - t(mu2[[i]])[1, ])^2) / (dim(X2)[1] * (dmax - d))
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * sigma2[[j]]^(-(dmax - d) / 2) * 
        exp(-(X[i, 1 : d, drop = FALSE] - t(mu1[[j]])) %*% 
              solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu1[[j]]) / 2 - 
              sum((X[i, (d + 1) : dmax] - t(mu2[[j]])[1, ])^2) / (2 * sigma2[[j]]))
    }
  }
  loglik <- -(dmax * n / 2) * log(2 * pi)
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  #print(loglik)
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu1[[j]] <- matrix(rep(0, d), d, 1)
      mu2[[j]] <- matrix(rep(0, dmax - d), (dmax - d), 1)
      for (i in 1 : n) {
        mu1[[j]] <- mu1[[j]] + t[i,j] * t(X[i, 1 : d, drop = FALSE])
        mu2[[j]] <- mu2[[j]] + t[i,j] * t(X[i, (d + 1) : dmax, drop = FALSE])
      }
      mu1[[j]] <- mu1[[j]] / sum(t[,j])
      mu2[[j]] <- mu2[[j]] / sum(t[,j])
      Sigma1[[j]] <- matrix(rep(0, d^2), d, d)
      sigma2[[j]] <- 0
      for (i in 1 : n) {
        Sigma1[[j]] <- Sigma1[[j]] + t[i,j] * (t(X[i, 1 : d, drop = FALSE]) - mu1[[j]]) %*% 
          (X[i, 1 : d, drop = FALSE] - t(mu1[[j]]))
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * sum((X[i, (d + 1) : dmax] - t(mu2[[j]])[1, ])^2)
      }
      Sigma1[[j]] <- Sigma1[[j]] / sum(t[,j])
      sigma2[[j]] <- sigma2[[j]] / ((dmax - d) * sum(t[,j]))
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * sigma2[[j]]^(-(dmax - d) / 2) * 
          exp(-(X[i, 1 : d, drop = FALSE] - t(mu1[[j]])) %*% 
                solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu1[[j]]) / 2 - 
                sum((X[i, (d + 1) : dmax] - t(mu2[[j]])[1, ])^2) / (2 * sigma2[[j]]))
      }
    }
    preLoglik <- loglik
    loglik <- - (dmax * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (dmax + (d + 1) * d / 2 + 1)
  bic <- 2 * loglik - numOfParameters * log(n)
  classification <- rep(0, n)
  for (i in 1 : n) {
    classification[i] <- which.max(f[i,])
  }
  return (list(bic = bic,
               classification = classification,
               d = d,
               k = k))
}

# calculate the bic for VVV model with number of clusters k, then do clustering
bicGmmVVV <- function(X, k, tao = 0) {
  # print("VVV")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu <- NULL
  Sigma <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), , drop = FALSE]
    head <- head + tao[i]
    mu[[i]] <- t(t(colMeans(X1)))
    Sigma[[i]] <- matrix(rep(0, d*d), d, d)
    for (l in 1 : dim(X1)[1]) {
      Sigma[[i]] <- Sigma[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*% 
        (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
    }
  }
  
  # Check if there is singularity
  for (j in 1 : k) {
    if (kappa(Sigma[[j]]) > 1e10) {
      print("Singular at initialization")
      return (list(bic = 0,
                   classification = init,
                   d = d,
                   k = k,
                   iter = 0,
                   logliks = NULL,
                   tao = tao,
                   mu = mu,
                   Sigma = Sigma))
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  fTemp <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      if (abs(myDet(Sigma[[j]])) < 1e-10) {
        f[i,j] <- tao[j] * (1e-10)^(-0.5) * 
          exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*%
                (1e10 * diag(d)) %*% (t(X[i, , drop = FALSE]) - mu[[j]]) / 2)
      } else {
        f[i,j] <- tao[j] * (myDet(Sigma[[j]])^(-0.5)) * 
          exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*%
                solve(Sigma[[j]]) %*% (t(X[i, , drop = FALSE]) - mu[[j]]) / 2)
      }
    }
  }
  loglik <- -(d * n / 2) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  
  flag <- 0
  while ((flag == 0) && (diff > 1e-5) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      if (sum(f[i,]) == 0) {
        for (j in 1 : k) {
          t[i,j] <- 1/k
        }
      } else {
        for (j in 1 : k) {
          t[i,j] <- f[i,j] / sum(f[i,])
          if (t[i,j] < 1e-10) {
            t[i,j] <- 1e-10
          }
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu[[j]] <- matrix(rep(0, d), d, 1)
      for (i in 1 : n) {
        mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, , drop = FALSE])
      }
      mu[[j]] <- mu[[j]] / sum(t[,j])
      Sigma[[j]] <- matrix(rep(0, d^2), d, d)
      for (i in 1 : n) {
        Sigma[[j]] <- Sigma[[j]] + t[i,j] * (t(X[i, , drop = FALSE]) - mu[[j]]) %*% 
          (X[i, , drop = FALSE] - t(mu[[j]]))
      }
      Sigma[[j]] <- Sigma[[j]] / sum(t[,j])
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        if (abs(myDet(Sigma[[j]])) < 1e-10) {
          flag <- 1
        } else {
          fTemp[i,j] <- tao[j] * (myDet(Sigma[[j]])^(-0.5)) * 
            exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*% 
                  solve(Sigma[[j]]) %*% (t(X[i, , drop = FALSE]) - mu[[j]]) / 2)
        }
      }
    }
    
    if (flag == 0) {
      f <- fTemp
    }
    
    preLoglik <- loglik
    loglik <- - (d * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2)
  bic <- 2 * loglik - numOfParameters * log(n)
  classification <- rep(0, n)
  for (i in 1 : n) {
    classification[i] <- which.max(f[i,])
  }
  return (list(bic = bic,
               loglik <- loglik,
               classification = classification))
}

# calculate the bic for diagonal model 1 with number of clusters k
# \mu_{1} = ... = \mu_{d} = 0; \sigma_{1}, ... , \sigma_{d} don't need to be equal
bicGmmDiag1 <- function(X, k, tao = 0) {
  # print("Diag1")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), , drop = FALSE]
    head <- head + tao[i]
    sigma2[[i]] <- rep(0, d)
    for (j in 1 : d) {
      for (l in 1 : dim(X1)[1]) {
        sigma2[[i]][j] <- sigma2[[i]][j] + X1[l, j]^2 / dim(X1)[1]
      }
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * prod(sigma2[[j]])^(-0.5) * 
        exp(-X[i, , drop = FALSE] %*% diag(1 / sigma2[[j]], d, d) %*%
              t(X[i, , drop = FALSE]) / 2)
    }
  }
  loglik <- -(d * n / 2) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      sigma2[[j]] <- rep(0, d)
      for (i in 1 : n) {
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * X[i, ]^2
      }
      sigma2[[j]] <- sigma2[[j]] / sum(t[,j])
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        f[i,j] <- tao[j] * prod(sigma2[[j]])^(-0.5) * 
          exp(-X[i, , drop = FALSE] %*% diag(1 / sigma2[[j]], d, d) %*%
                t(X[i, , drop = FALSE]) / 2)
      }
    }
    preLoglik <- loglik
    loglik <- - (d * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * d
  bic <- 2 * loglik - numOfParameters * log(n)
  return (list(bic = bic))
}

# calculate the bic for diagonal model 2 with number of clusters k
# \mu_{1} = ... = \mu_{d} = 0; \sigma_{1} = ... = \sigma_{d}
bicGmmDiag2 <- function(X, k, tao = 0) {
  # print("Diag2")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), , drop = FALSE]
    head <- head + tao[i]
    sigma2[[i]] <- 0
    for (l in 1 : dim(X1)[1]) {
      sigma2[[i]] <- sigma2[[i]] + sum(X1[l, ]^2) / (dim(X1)[1] * d)
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  fTemp <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * sigma2[[j]]^(-d / 2) * exp(-sum(X[i, ]^2) / (2 * sigma2[[j]]))
    }
  }
  loglik <- -(d * n / 2) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  
  flag <- 0
  while ((flag == 0) && (diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      sigma2[[j]] <- 0
      for (i in 1 : n) {
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * sum(X[i, ]^2)
      }
      sigma2[[j]] <- sigma2[[j]] / (d * sum(t[,j]))
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        if (abs(sigma2) < 1e-10) {
          flag <- 1
        }
        fTemp[i,j] <- tao[j] * sigma2[[j]]^(-d / 2) * exp(-sum(X[i, ]^2) / (2 * sigma2[[j]]))
      }
    }
    if (flag == 0) {
      f <- fTemp
    }
    preLoglik <- loglik
    loglik <- - (d * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k
  bic <- 2 * loglik - numOfParameters * log(n)
  return (list(bic = bic))
}

# calculate the bic for diagonal model 3 with number of clusters k
# \mu_{1}, ... , \mu_{d} have no contraints; \sigma_{1}, ... , \sigma_{d} don't need to be equal
bicGmmDiag3 <- function(X, k, tao = 0) {
  # print("Diag3")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu <- NULL
  sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), , drop = FALSE]
    head <- head + tao[i]
    mu[[i]] <- t(t(colMeans(X1)))
    sigma2[[i]] <- rep(0, d)
    for (j in 1 : d) {
      for (l in 1 : dim(X1)[1]) {
        sigma2[[i]][j] <- sigma2[[i]][j] + (X1[l, j] - mu[[i]][j])^2 / dim(X1)[1]
      }
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * prod(sigma2[[j]])^(-0.5) * 
        exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*% diag(1 / sigma2[[j]], d, d) %*%
              (t(X[i, , drop = FALSE]) - mu[[j]])/ 2)
    }
  }
  loglik <- -(d * n / 2) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu[[j]] <- matrix(rep(0, d), d, 1)
      for (i in 1 : n) {
        mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, , drop = FALSE])
      }
      mu[[j]] <- mu[[j]] / sum(t[,j])
      sigma2[[j]] <- rep(0, d)
      for (i in 1 : n) {
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * (X[i, ] - t(mu[[j]])[1, ])^2
      }
      sigma2[[j]] <- sigma2[[j]] / sum(t[,j])
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        f[i,j] <- tao[j] * prod(sigma2[[j]])^(-0.5) * 
          exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*% diag(1 / sigma2[[j]], d, d) %*%
                (t(X[i, , drop = FALSE]) - mu[[j]])/ 2)
      }
    }
    preLoglik <- loglik
    loglik <- - (d * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (d + d)
  bic <- 2 * loglik - numOfParameters * log(n)
  return (list(bic = bic))
}

# calculate the bic for diagonal model 4 with number of clusters k
# \mu_{1}, ... , \mu_{d} have no contraints; \sigma_{1} = ... = \sigma_{d}
bicGmmDiag4 <- function(X, k, tao = 0) {
  # print("Diag4")
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu <- NULL
  sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), , drop = FALSE]
    head <- head + tao[i]
    mu[[i]] <- t(t(colMeans(X1)))
    sigma2[[i]] <- 0
    for (l in 1 : dim(X1)[1]) {
      sigma2[[i]] <- sigma2[[i]] + sum((X1[l, ] - t(mu[[i]])[1, ])^2) / (dim(X1)[1] * d)
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      f[i,j] <- tao[j] * sigma2[[j]]^(-d / 2) * 
        exp(-sum((X[i, ] - t(mu[[j]])[1, ])^2) / (2 * sigma2[[j]]))
    }
  }
  loglik <- -(d * n / 2) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu[[j]] <- matrix(rep(0, d), d, 1)
      for (i in 1 : n) {
        mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, , drop = FALSE])
      }
      mu[[j]] <- mu[[j]] / sum(t[,j])
      sigma2[[j]] <- 0
      for (i in 1 : n) {
        sigma2[[j]] <- sigma2[[j]] + t[i,j] * sum((X[i, ] - t(mu[[j]])[1, ])^2)
      }
      sigma2[[j]] <- sigma2[[j]] / (d * sum(t[,j]))
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        f[i,j] <- tao[j] * sigma2[[j]]^(-d / 2) * 
          exp(-sum((X[i, ] - t(mu[[j]])[1, ])^2) / (2 * sigma2[[j]]))
      }
    }
    preLoglik <- loglik
    loglik <- - (d * n / 2) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (d + 1)
  bic <- 2 * loglik - numOfParameters * log(n)
  return (list(bic = bic))
}

# calulate the BIC of the blockwise model 1 with 2*2 block after d dimension, and then do clustering
# \mu_{d+1} = ... = \mu_{dmax} = 0; 2*2 \Sigma's don't need to be equal
bicGmmPair1 <- function(X, d, k, tao = 0) {
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  dmax <- dim(X)[2] / 2
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu <- NULL
  Sigma1 <- NULL
  Sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), 1 : (2*d), drop = FALSE]
    X2 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), (2*d + 1) : (2*dmax), drop = FALSE]
    head <- head + tao[i]
    mu[[i]] <- t(t(colMeans(X1)))
    Sigma1[[i]] <- matrix(rep(0, (2*d)*(2*d)), 2*d, 2*d)
    for (l in 1 : dim(X1)[1]) {
      Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*%
        (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
    }
    Sigma2[[i]] <- NULL
    for (j in 1 : (dmax - d)) {
      Sigma2[[i]][[j]] <- matrix(rep(0,4),2,2)
      for (l in 1 : dim(X2)[1]) {
        Sigma2[[i]][[j]] <- Sigma2[[i]][[j]] + t(X2[l, (2*j-1) : (2*j), drop = FALSE]) %*%
          X2[l, (2*j-1) : (2*j), drop = FALSE] / dim(X2)[1]
      }
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      prodTerm <- 1
      sumTerm <- 0
      for (l in 1 : (dmax - d)) {
        prodTerm <- prodTerm * myDet(Sigma2[[j]][[l]])^(-0.5)
        if (sum(abs(Sigma2[[j]][[l]])) == 0) {
          sumTerm <- sumTerm + 1e10
        } else {
          sumTerm <- sumTerm + X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE] %*% 
            solve(Sigma2[[j]][[l]]) %*% t(X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE])
        }
      }
      f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * prodTerm * 
        exp(-(X[i, 1 : (2*d), drop = FALSE] - t(mu[[j]])) %*% 
              solve(Sigma1[[j]]) %*% (t(X[i, 1 : (2*d), drop = FALSE]) - mu[[j]]) / 2 - 
              sumTerm / 2)
    }
  }
  loglik <- -(dmax * n) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  #print(loglik)
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu[[j]] <- matrix(rep(0, (2*d)), 2*d, 1)
      for (i in 1 : n) {
        mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, 1 : (2*d), drop = FALSE])
      }
      mu[[j]] <- mu[[j]] / sum(t[,j])
      Sigma1[[j]] <- matrix(rep(0, (2*d)*(2*d)), 2*d, 2*d)
      for (i in 1 : n) {
        Sigma1[[j]] <- Sigma1[[j]] + t[i,j] * (t(X[i, 1 : (2*d), drop = FALSE]) - mu[[j]]) %*% 
          (X[i, 1 : (2*d), drop = FALSE] - t(mu[[j]]))
      }
      for (l in 1 : (dmax - d)) {
        Sigma2[[j]][[l]] <- matrix(rep(0,4),2,2)
        for (i in 1 : n) {
          Sigma2[[j]][[l]] <- Sigma2[[j]][[l]] + t[i,j] * 
            t(X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE]) %*% 
            (X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE])
        }
        Sigma2[[j]][[l]] <- Sigma2[[j]][[l]] / sum(t[,j])
      }
      Sigma1[[j]] <- Sigma1[[j]] / sum(t[,j])
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        prodTerm <- 1
        sumTerm <- 0
        for (l in 1 : (dmax - d)) {
          prodTerm <- prodTerm * myDet(Sigma2[[j]][[l]])^(-0.5)
          sumTerm <- sumTerm + X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE] %*% 
            solve(Sigma2[[j]][[l]]) %*% t(X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE])
        }
        f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * prodTerm * 
          exp(-(X[i, 1 : (2*d), drop = FALSE] - t(mu[[j]])) %*% 
                solve(Sigma1[[j]]) %*% (t(X[i, 1 : (2*d), drop = FALSE]) - mu[[j]]) / 2 - 
                sumTerm / 2)
      }
    }
    preLoglik <- loglik
    loglik <- - (dmax * n) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (2*d + (2*d + 1) * (2*d) / 2 + 3*(dmax - d))
  bic <- 2 * loglik - numOfParameters * log(n)
  classification <- rep(0, n)
  for (i in 1 : n) {
    classification[i] <- which.max(f[i,])
  }
  return (list(bic = bic,
               classification = classification,
               d = d,
               k = k))
}

# calulate the BIC of the blockwise model 1 with 2*2 block after d dimension, and then do clustering
# \mu_{d+1} = ... = \mu_{dmax} = 0; 2*2 \Sigma's must be equal
bicGmmPair2 <- function(X, d, k, tao = 0) {
  diff <- 1
  iter <- 0
  n <- dim(X)[1]
  dmax <- dim(X)[2] / 2
  
  # initialization of parameters
  if ((sum(tao) != 1) || (length(tao) != k)) {
    tao <- rep(1 / k, k)
  }
  mu <- NULL
  Sigma1 <- NULL
  Sigma2 <- NULL
  head <- 0
  for (i in 1 : k) {
    X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), 1 : (2*d), drop = FALSE]
    X2 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), (2*d + 1) : (2*dmax), drop = FALSE]
    head <- head + tao[i]
    mu[[i]] <- t(t(colMeans(X1)))
    Sigma1[[i]] <- matrix(rep(0, (2*d)*(2*d)), 2*d, 2*d)
    for (l in 1 : dim(X1)[1]) {
      Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*% 
        (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
    }
    Sigma2[[i]] <- matrix(rep(0,4),2,2)
    for (j in 1 : (dmax - d)) {
      for (l in 1 : dim(X2)[1]) {
        Sigma2[[i]] <- Sigma2[[i]] + t(X2[l, (2*j-1) : (2*j), drop = FALSE]) %*% 
          X2[l, (2*j-1) : (2*j), drop = FALSE] / (dim(X2)[1]*(dmax - d))
      }
    }
  }
  
  # initialization of f and t
  f <- matrix(rep(0, n * k), n, k)
  t <- matrix(rep(0, n * k), n, k)
  for (i in 1 : n) {
    for (j in 1 : k) {
      sumTerm <- 0
      for (l in 1 : (dmax - d)) {
        if (sum(abs(Sigma2[[j]])) == 0) {
          sumTerm <- sumTerm + 1e10
        } else {
          sumTerm <- sumTerm + X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE] %*% 
            solve(Sigma2[[j]]) %*% t(X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE])
        }
      }
      f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * myDet(Sigma2[[j]])^(-(dmax-d)/2) * 
        exp(-(X[i, 1 : (2*d), drop = FALSE] - t(mu[[j]])) %*% 
              solve(Sigma1[[j]]) %*% (t(X[i, 1 : (2*d), drop = FALSE]) - mu[[j]]) / 2 - 
              sumTerm / 2)
    }
  }
  loglik <- -(dmax * n) * log(2 * pi) 
  for (i in 1 : n) {
    loglik <- loglik + log(sum(f[i,]))
  }
  #print(loglik)
  
  while ((diff > 1e-2) && (iter < 100)) {
    # E step
    for (i in 1 : n) {
      for (j in 1 : k) {
        t[i,j] <- f[i,j] / sum(f[i,])
        if (t[i,j] < 1e-10) {
          t[i,j] <- 1e-10
        }
      }
    }
    # M step
    for (j in 1 : k) {
      tao[j] <- sum(t[,j]) / n
      mu[[j]] <- matrix(rep(0, (2*d)), 2*d, 1)
      for (i in 1 : n) {
        mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, 1 : (2*d), drop = FALSE])
      }
      mu[[j]] <- mu[[j]] / sum(t[,j])
      Sigma1[[j]] <- matrix(rep(0, (2*d)*(2*d)), 2*d, 2*d)
      Sigma2[[j]] <- matrix(rep(0,4),2,2)
      for (i in 1 : n) {
        Sigma1[[j]] <- Sigma1[[j]] + t[i,j] * (t(X[i, 1 : (2*d), drop = FALSE]) - mu[[j]]) %*% 
          (X[i, 1 : (2*d), drop = FALSE] - t(mu[[j]]))
      }
      for (l in 1 : (dmax - d)) {
        for (i in 1 : n) {
          Sigma2[[j]] <- Sigma2[[j]] + t[i,j] * 
            t(X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE]) %*% 
            (X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE])
        }
      }
      
      Sigma1[[j]] <- Sigma1[[j]] / sum(t[,j])
      Sigma2[[j]] <- Sigma2[[j]] / ((dmax-d)*sum(t[,j]))
    }
    # check convergence
    for (i in 1 : n) {
      for (j in 1 : k) {
        sumTerm <- 0
        for (l in 1 : (dmax - d)) {
          sumTerm <- sumTerm + X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE] %*% 
            solve(Sigma2[[j]]) %*% t(X[i, (2*d+2*l-1) : (2*d+2*l), drop = FALSE])
        }
        f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * myDet(Sigma2[[j]])^(-(dmax-d)/2) * 
          exp(-(X[i, 1 : (2*d), drop = FALSE] - t(mu[[j]])) %*% 
                solve(Sigma1[[j]]) %*% (t(X[i, 1 : (2*d), drop = FALSE]) - mu[[j]]) / 2 - 
                sumTerm / 2)
      }
    }
    preLoglik <- loglik
    loglik <- - (dmax * n) * log(2 * pi)
    for (i in 1 : n) {
      loglik <- loglik + log(sum(f[i,]))
    }
    #print(loglik)
    diff <- abs(loglik - preLoglik)
    iter <- iter + 1
  }
  # calculate BIC
  numOfParameters <- (k - 1) + k * (2*d + (2*d + 1) * (2*d) / 2 + 3)
  bic <- 2 * loglik - numOfParameters * log(n)
  classification <- rep(0, n)
  for (i in 1 : n) {
    classification[i] <- which.max(f[i,])
  }
  return (list(bic = bic,
               classification = classification,
               d = d,
               k = k))
}

myDet <- function(X) {
  if (is.matrix(X)) {
    result <- det(X)
  } else {
    result <- X[1]
  }
  return (result)
}