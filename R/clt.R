library("irlba")
library("mclust")
library("cluster")
library("MASS")
library("reshape")
library("plyr")
library("ggplot2")

## Sample an undirected graph on n vertices
## Input: P is n times n matrix giving the parameters of the Bernoulli r.v.
rg.sample <- function(P){
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(A)
}

rg.sample.pois <- function(P){
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- rpois(n*(n-1)/2, lambda = P[col(P) > row(P)])
  U <-  U + t(U)
  return(U)
}

rg.SBM <- function(n, B, rho,condition=FALSE){
  if(!condition){
    tau <- sample(c(1:length(rho)), n, replace = TRUE, prob = rho)
  }
  else{
    tau <- unlist(lapply(1:2,function(k) rep(k, rho[k]*n)))
  }
  P <- B[tau,tau]
  return(list(adjacency=rg.sample(P),tau=tau))
}

stfp <- function(A, dim, method = "svd"){

    if(method == "svd"){
        S <- irlba(A, nu = dim, nv = dim)
        V <- S$v[,1:dim]
        D <- S$d[1:dim]
    } else {
        S <- eigen(A)
        V <- S$vectors[,1:dim]
        D <- S$values[1:dim]
    }

    if(dim == 1){
        Xhat <- V*sqrt(D)
    }
    else{
      Xhat <- V %*% diag(sqrt(D))
    }

    return(Xhat)
}

procrustes <- function(A,B){
    mm <- svd(t(A) %*% B)
    W <- mm$u %*% t(mm$v)
    return(W)
}

sd.multivariate <- function(X){
    n <- nrow(X)
    Xbar <- outer(rep(1,n), colMeans(X))
    Xtilde <- X - Xbar
    return(1/(n-1)*(t(Xtilde)%*%Xtilde))
}

block.var <- function(X,rho){
    n <-  nrow(X)
    var.list <- list()

    Delta <-  matrix(0, nrow = ncol(X), ncol = ncol(X))
    for( i in 1:n){
        Delta <- Delta + outer(X[i,],X[i,])*rho[i]
    }

    for(i in 1:n){
        tmp1 <- X[i,]%*%t(X)
        tmp2 <- tmp1 - tmp1^2
        B <- matrix(0, nrow = ncol(X), ncol = ncol(X))
        for(j in 1:n){
            B <-  B + outer(X[j,], X[j,])*tmp2[j]*rho[j]
        }
        var.list[[i]] <- solve(Delta)%*%B%*%solve(Delta)
    }
    return(var.list)
}

bayes.mvn <- function(n, mu1, mu2, Sigma1, Sigma2, pi1, pi2, nmc){

    tau <- sample(1:2, nmc, replace = TRUE, prob = c(pi1,pi2))
    m1 <- sum(tau == 1)
    m2 <-  sum(tau == 2)

    X1 <- mvrnorm(m1, sqrt(n)*mu1, Sigma1)
    X2 <-  mvrnorm(m2, sqrt(n)*mu2, Sigma2)
    labels <- c(rep(1,m1),rep(-1,m2))

    x <- rbind(X1,X2)

    Sigma1.inv <- solve(Sigma1)
    Sigma2.inv <- solve(Sigma2)

    c1 <- log(pi1) - 1/2*log(det(Sigma1))
    c2 <- log(pi2) - 1/2*log(det(Sigma2))

    xbar1 <- x - sqrt(n)*outer(rep(1,nrow(x)),mu1)
    xbar2 <- x - sqrt(n)*outer(rep(1,nrow(x)),mu2)

    density1 <- c1 - rowSums((1/2*xbar1%*%Sigma1.inv)*xbar1)
    density2 <- c2 - rowSums((1/2*xbar2%*%Sigma2.inv)*xbar2)

    labels.hat <- sign(density1 - density2)
    error <- sum(labels.hat != labels)/nmc
    return(error)
}


experiment1 <- function(n = 1000){
    #n <- 8000
    rho <- c(0.6,0.4)
    B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)

    A <- rg.SBM(n, B,rho)
    tau <- A$tau
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))
    X <- xx[tau,]

    Xhat <- stfp(A$adjacency,dim = 2)
    W <- procrustes(Xhat,X)
    residue <- sqrt(n)*(Xhat%*%W - X)
    print(norm(Xhat%*%W - X,"F")^2)

    sigs <- block.var(xx,rho)

    print(sigs)

   #  mean.residue1 <- colMeans(residue[which(tau == 1),])
   sd.residue1 <- sd.multivariate(residue[which(tau == 1),])
   #  mean.residue2 <- colMeans(residue[which(tau == 2),])
   sd.residue2 <- sd.multivariate(residue[which(tau == 2),])
   print(sd.residue1)
   print(sd.residue2)

   # ## pdf(paste("clusplot",n,".pdf",sep=''),useDingbats=FALSE)
   #  aa <- Mclust(sqrt(n)*Xhat,2,modelNames = c("VVV"))
   #  error.rate.gmm <- min(sum(abs(aa$classification - tau)),
   #                    sum(abs(3 - aa$classification - tau)))/n
    df <- data.frame(x = (Xhat%*%W)[,1], y = (Xhat%*%W)[,2],
                     block = as.character(tau),
                     tau = tau)
##    qplot(data = dat, x = X1, y = X2, colour = labels) + stat_ellipse()

    library(ellipse)
    df_ell <- data.frame()

    for(g in unique(df$tau)){
        sig <- sigs[[g]]
        df_ell <- rbind(df_ell, cbind(as.data.frame(
          with(df[df$tau==g,],
               ellipse(sig[1,2]/sqrt(sig[1,1]*sig[2,2]),
               scale=c(sqrt(sig[1,1])/sqrt(n),sqrt(sig[2,2])/sqrt(n)),
               centre=c(xx[g,1],xx[g,2]))),group=g)))

    }

    library(ggplot2)
    ggplot(data=df, aes(x=x, y=y,color=block)) +
      geom_point(size=2) +  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
      geom_path(data=df_ell[1:(nrow(df_ell)/2),], aes(x=x, y=y), size=1,linetype = 2, colour=1) +
          geom_path(data = df_ell[-c(1:(nrow(df_ell)/2)),], aes(x = x, y = y), size = 1, linetype = 2, colour = 1) + theme(legend.position = "none")
    ggsave(paste("clusplot",n,".pdf",sep=''),width=6,height=4)


    # clusplot(sqrt(n)*Xhat%*%W,aa$classification, color=TRUE, shade = TRUE,
    #          span = FALSE, sub = "",
    #          main = paste("Gaussian mixture clustering of the residue for n = ",
    #            n, sep = ''), xlab = "", ylab = "")
    # dev.off()

    # bb <- kmeans(sqrt(n)*Xhat,2, iter.max = 50)
    # clusplot(sqrt(n)*Xhat%*%W,bb$cluster,color = TRUE, shade = TRUE, sub = "", main = paste("Kmeans clustering of the residue for n = ", n, sep = ''))

    # error.rate.kmeans <- min(sum(abs(bb$cluster - tau)),
    #                   sum(abs(3 - bb$cluster - tau)))/n

    ## return(df)

}

# B2 <- rbind(c(0.2,0.1),c(0.1,0.25))
# B3 <- rbind(c(0.2,0.1,0.08),c(0.1,0.25,0.05),
#             c(0.08,0.05,0.4))
# pi2 <- c(0.5,0.5)
# pi3 <- c(0.4,0.4,0.2)

experiment1_YP <- function(n = 1000){
  #n <- 8000
  set.seed(12345)
#  rho <- c(0.5,0.5)
#  B <- matrix(c(0.2,0.1,0.1,0.25), nrow = 2, ncol = 2)
  rho <- c(0.4,0.4,0.2)
  B <- rbind(c(0.2,0.1,0.08),c(0.1,0.25,0.05), c(0.08,0.05,0.4))
  K <- length(rho)

  A <- rg.SBM(n, B,rho)
  tau <- A$tau
  x <- eigen(B)
  xx <- x$vectors %*% diag(sqrt(x$values))
  X <- xx[tau,]

  Xhat <- stfp(A$adjacency, dim=K)
  W <- procrustes(Xhat,X)
  residue <- sqrt(n)*(Xhat%*%W - X)
  print(norm(Xhat%*%W - X,"F")^2)

  sigs <- block.var(xx,rho)
  print(sigs)

  # #  mean.residue1 <- colMeans(residue[which(tau == 1),])
  # sd.residue1 <- sd.multivariate(residue[which(tau == 1),])
  # #  mean.residue2 <- colMeans(residue[which(tau == 2),])
  # sd.residue2 <- sd.multivariate(residue[which(tau == 2),])
  # print(sd.residue1)
  # print(sd.residue2)

  df <- data.frame(x = (Xhat%*%W)[,1], y = (Xhat%*%W)[,2],
                   block = as.character(tau),
                   tau = tau)

  library(ellipse)
  df_ell <- data.frame()
  for(g in unique(df$tau)){
    sig <- sigs[[g]]
    df_ell <- rbind(df_ell, cbind(as.data.frame(
      with(df[df$tau==g,],
           ellipse(sig[1,2]/sqrt(sig[1,1]*sig[2,2]),
                   scale=c(sqrt(sig[1,1])/sqrt(n),sqrt(sig[2,2])/sqrt(n)),
                   centre=c(xx[g,1],xx[g,2]))),group=g)))

  }

  df_ell <- cbind(df_ell, block=as.character(rep(1:K, each=nrow(df_ell)/K)))
  df.xx <- as.data.frame(xx[,1:2]); names(df.xx) <- c("x","y")
  df.xx <- cbind(df.xx, block=as.character(c(1:K)))

  library(ggplot2)
  ggplot(data=df, aes(x=x, y=y, color=block)) +
    geom_point(size=2, alpha=0.5) +
    stat_ellipse(type="norm", size=1, linetype = "dashed", aes(color=block)) +
    geom_point(data=df.xx, aes(x=x,y=y), color="black", size=3) +
    geom_path(data=df_ell, aes(x=x, y=y, group=block), size=1, linetype = 2, colour=1) +
#    geom_path(data = df_ell[-c(1:(nrow(df_ell)/K)),], aes(x = x, y = y), size = 1, linetype = 2, colour = 1) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  #  theme_light() +
    theme(legend.position = "none")
  ggsave(paste0("Fig1-n",n,"-K",K,".pdf"),width=5,height=5)
}

experiment2 <- function(){
    set.seed(12345)
    nseq <- seq(from = 1000, to = 4000, by = 250)
    rho <- c(0.5,0.5)
    B <- matrix(c(0.62,0.32,0.32,0.62), nrow = 2, ncol = 2)
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))

    error.vec.gmm <- numeric(length(nseq))
    error.vec.kmeans <-  numeric(length(nseq))
    error.vec.bayes <- numeric(length(nseq))

    tmp <- block.var(xx, rho)
    Sigma1 <- tmp[[1]]
    Sigma2 <- tmp[[2]]
    nmc <- 100

    for(i in 1:length(nseq)){
        error.i.gmm <- numeric(nmc)
        error.i.kmeans <- numeric(nmc)
        for(j in 1:nmc){
            ## A <- rg.SBM(nseq[i], B, rho)
            ##tau <- rep(c(1:nrow(B)),nseq[i]*rho)
            ## tau <- sample(1:nrow(B), nseq[i], replace = TRUE, prob = rho)
            ## P <- B[tau,tau]
            ## A <- rg.sample(P)
            ## X <- xx[tau,]
            rho.conditional <- rmultinom(1, nseq[i], rho)
            g <- sbm.game(nseq[i], B, rho.conditional)
            tau <- rep(1:length(rho), rho.conditional)

##            diag(A) <- rowSums(A)/sqrt(sum(rowSums(A)))
##            diag(A) <- rowSums(A)/(nseq[i] - 1)

#            Xhat <- (adjacency.spectral.embedding(g, 4)$X)[,1:2]
            Xhat <- (embed_adjacency_matrix(g, 4)$X)[,1:2]
            aa <- Mclust(Xhat,2,modelNames = c("VVV"))
            tmp <- sum(aa$classification != tau)/nseq[i]
            error.i.gmm[j] <- min(tmp, 1 - tmp)

            bb <- kmeans(Xhat,2, iter.max = 50)
            tmp <- sum(bb$cluster != tau)/nseq[i]
            error.i.kmeans[j] <- min(tmp, 1 - tmp)
        }
        error.vec.gmm[i] <- mean(error.i.gmm)
        error.vec.kmeans[i] <- mean(error.i.kmeans)
        error.vec.bayes[i] <- bayes.mvn(nseq[i], xx[1,], xx[2,], Sigma1, Sigma2, rho[1], rho[2], 100000)
   }
        error.vec.log <- log(nseq)/nseq
        error.vec.log <- error.vec.kmeans[1]/error.vec.log[1]*error.vec.log

        dat.new <- melt(data.frame(nseq, gmm = log10(error.vec.gmm), kmeans = log10(error.vec.kmeans), bayes = log10(error.vec.bayes), log = log10(error.vec.log)), id = "nseq")
        ggplot(dat.new, aes(x = nseq, y = value, colour = variable)) + geom_line() + xlab("n") +
          ylab("classification error") + labs(title = element_blank()) + theme_light() + theme(axis.text = element_text(size = 15)) +
            ## theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
                ## theme(legend.key = element_rect(fill = 'white', colour = 'white')) + scale_color_discrete(breaks=c("gmm", "kmeans", "bayes", "log"), labels = c("GMM", "K-Means", "Bayes", "STFP"))
       ggsave("gmm_kmeans_bayes.pdf", width=6,height=4)

}

# if sum(rho) == 1, randomly allocated block memberships
# if sum(rho) == n, conditionally allocated block memberships
sbm = function(n,K,B,rho,myseed=F){
if(myseed) set.seed(myseed)
if(sum(rho) == n) blockn = rho
if(sum(rho) != n) blockn = tabulate(sample(K,n,replace=T,prob=rho),K)
tau=NULL; for(k in 1:K) tau = c(tau,rep(k,blockn[k]))
A = matrix(0,nrow=n,ncol=n)
for(i in 1:(n-1)) for(j in (i+1):n) A[i,j] = A[j,i] = rbinom(1,1,B[tau[i],tau[j]])
return(list(A,tau))
}

experiment2.dan <- function(nmc=100, nseq=seq(from = 1000, to = 4000, by = 250)){
  rho <- c(0.6,0.4)
  B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)
  x <- eigen(B)
  xx <- x$vectors %*% diag(sqrt(x$values))

  tmp <- block.var(xx, rho)
  Sigma1 <- tmp[[1]]
  Sigma2 <- tmp[[2]]


  all <- data.frame(n=numeric(0),error=numeric(0),name=character(0))
  error.max = 0
  for(n in nseq){

    for(j in 1:nmc){
      A <- rg.SBM(n, B, rho)
      tau <- A$tau
      X <- xx[tau,]

      Xhat <- stfp(A$adjacency,dim = 2, method="eigen")
      aa <- Mclust(sqrt(n)*Xhat,2,modelNames = c("VVV"))
      error.gmm <- min(sum(abs(aa$classification - tau)),
                            sum(abs(3 - aa$classification - tau)))/n

      bb <- kmeans(sqrt(n)*Xhat,2, iter.max = 50)
      error.kmeans <- min(sum(abs(bb$cluster - tau)),
                               sum(abs(3 - bb$cluster - tau)))/n

      error.max = max(error.max,error.kmeans)

      all <- rbind(all,data.frame(n=c(n,n), err = c(error.gmm,error.kmeans),name=c("gmm","kmeans")))

    }
    error.bayes <- bayes.mvn(n, xx[1,], xx[2,], Sigma1, Sigma2, rho[1], rho[2], 1000000)
    error.log <- error.max/(log(nseq[1])/nseq[1])*log(n)/n

    all <- rbind(all, data.frame(n=c(n,n),
                                 err = c(error.bayes,error.log),
                                 name=c("oracle bayes","log bound")))
  }

  all.stat <- ddply(all,.(n,name),
    function(nn) data.frame(mean=mean(nn$err), se=sd(nn$err)/sqrt(nmc)))
  ggplot(all.stat)+aes(x=n,y=mean,ymin=mean-2*se,ymax=mean+2*se,color=name)+
    geom_line()+
    geom_ribbon(aes(linetype=NA),alpha=.2)+
    scale_y_log10()
  ggsave("gmm_kmeans_bayes.pdf", width=6,height=4)

  all
}

experiment3.dan <- function(nmc=100, nseq=seq(from = 1000, to = 4000, by = 250)){
  start.time <- Sys.time()
  rho <- c(0.6,0.4)
  B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)
  x <- eigen(B)
  xx <- x$vectors %*% diag(sqrt(x$values))

  all <- data.frame(n=numeric(0),error=numeric(0),name=character(0))
  error.max = 0
  for(n in nseq){
    cat('n=',n,'\n')
    for(j in 1:nmc){
      A <- rg.SBM(n, B, rho,condition=TRUE)
      tau <- A$tau
      X <- xx[tau,]

      Xhat <- stfp(A$adjacency,dim = 2, method="eigen")

      W <- procrustes(Xhat,X)
      sq.error <- norm(Xhat%*%W - X,"f")^2
      all = rbind(all,data.frame(n=n,sq.error=sq.error,name="Observed"))
      all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(Xhat)),name="Predicted"))
      cat('\rMC = ',j,' ',format(Sys.time()-start.time))
    }
    all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(X)),name="Expected"))
    cat('\n')

    p<- ggplot(all)+aes(factor(n),sq.error)+geom_violin()
    ggsave("L2_Error_Violin.pdf",width=6,height=4)
  }

  all
}

expectation.apv <- function(X){
  n <- dim(X)[1]
  SP = diag(t(X) %*% X)
  V = X %*% diag(SP**(-1/2))
  P = X %*% t(X)
  VP =  P*(1-P)-diag(diag(P*(1-P)))
  E <- sum(VP %*% diag(X %*% t(X)))
  for(i in 1:n){
    E <- E+norm(as.matrix(X[i,]))**6 # * norm(as.matrix(V[i,]))**2
#     for(j in 1:n){
#       if(j != i){
#         E <- E+VP[i,j]*norm(as.matrix(X[j,]))**2
#       }
#     }
  }
  E / min(SP**2)
}

binomial.test <- function(p,q,r,s,n1,n2,nmc){

    X1 <- rbinom(0.6*nmc,n1,p)
    X2 <- rbinom(0.6*nmc,n2,q)

    Y1 <- rbinom(0.4*nmc,n1,r)
    Y2 <- rbinom(0.4*nmc,n2,s)

    Z1 <- c(X1,Y1)
    Z2 <- c(X2,Y2)
    labels <- c(rep(1,0.6*nmc),rep(-1,0.4*nmc))

    tt1 <- Z1*log(p) + (n1-Z1)*log(1-p) + Z2*log(q) + (n2-Z2)*log(1 - q)
    tt2 <- Z1*log(r) + (n1-Z1)*log(1-r) + Z2*log(s) + (n2-Z2)*log(1 - s)

    tt <- tt1 - tt2

    tt[abs(tt) <= 1e-14] <- -1
    yy <- sign(tt)

    error.bayes <- 1 - sum(sign(yy) == labels)/(nmc)
    return(error.bayes)
}

test1 <- function(n,nmc,p){
    tvec <- numeric(nmc)
    u <- rep(c(1:2),n/2)
    for(i in 1:nmc){
        P <- matrix(p, nrow = n, ncol = n)
        A <- rg.sample(P)
        Xhat <- stfp(A, dim = 1)
        tmp1 <- sum((Xhat - sqrt(p))^2)
        tmp2 <- sum((-Xhat - sqrt(p))^2)
        if(tmp1 < tmp2){
            mult <- 1
        } else {
            mult <- -1
        }

        Xhat <- Xhat*mult
        tvec[i] <- sum(u*(Xhat - sqrt(p)))
    }
    return(tvec)
}

example1 <- function(nmc,k){
    xxx <- numeric(nmc)

    for(i in 1:nmc){
        X <- matrix(runif(3*k,max=1/sqrt(3)),ncol = 3)
        rho <- rep(1,k)/k
        sigma.list <- block.var(X, rho)

        p1 <- sum(X[1,]^2)

        xxx[i] <- 2*(t(X[1,]) %*% sigma.list[[1]] %*% X[1,])/(p1*(1 - p1))
    }
    return(xxx)
}

p2q2 <- function(rho, p, q){
  delta <- rho*p^2 + (1 - rho)*q^2

  # part1 <- rho^2*p^2*(rho^2*p^6*(1 - p^2) + 2*rho*(1 - rho)*p^3*q^3*(1 - p*q) +
  #                     (1 - rho)^2*q^6*(1 - q^2))
  #
  # part2 <- 4*rho^2*p^4*(1 - p^2) + rho*(1 - rho)*p*q*(1 - p*q)*(p^2 + q^2)
  #
  # part3 <- 4*rho^3*p^6*(1 - p^2) + 2*rho^2*(1 - rho)*p^3*q^3*(1 - p*q) +
  #                                            2*rho^2*(1 - rho)*p^4*q^2*(1 - p*q)

  part1 <- rho^2*p^4*(1 - p^2)*(4 + rho^2*p^4/delta^2 - 4*rho*p^2/delta)/(2*delta^2)
  part2 <- 1/2*rho^2*(1 - rho)^2*p^2*q^6*(1 - q^2)/delta^4
  part3 <- rho*(1 - rho)*p*q*(1 - p*q)/delta*(q^2/delta + rho^2*p^4*q^2/delta^3 - 2*p^2*q^2*rho/delta^2)

  # part1 <- part1/(2*delta^4)
  # part2 <- part2/(2*delta^2)
  # part3 <- part3/(2*delta^3)

  pq.var <- part1 + part2 + part3


  naive.numerator <- rho*p^4*(1 - p^2) + (1 - rho)*p*q^3*(1 - p*q)

  pq.var2 <- rho*naive.numerator/delta^2

  a1 <- (1 - rho*p^2/(2*delta))
  a2 <- rho*p*q/(2*delta)

  part1 <- 2*rho^2*p^2*(1 - p^2)/delta*a1^2*p^2/delta
  part2a <- 0.5*rho*(1 - rho)*p*q*(1 - p*q)/delta*(a1^2*q^2/delta + a2^2*p^2/delta - 2*a1*a2*p*q/delta)
  part2b <- 0.5*rho*(1 - rho)*p*q*(1 - p*q)/delta*(a1^2*p^2/delta + a2^2*q^2/delta - 2*a1*a2*p*q/delta)
  part3 <- 2*(1 - rho)^2*q^2*(1 - q^2)/delta*a2^2*q^2/delta

  pq.var3 <- part1 + 2*part2a  + part3

  return(c(pq.var, pq.var2,pq.var3))
}

# nmc <- 1000
# phat <- numeric(nmc)
# qhat <- numeric(nmc)
# dphat <- numeric(nmc)
# dqhat <- numeric(nmc)
# p <- 0.42
# q <- 0.54
# n <- 2000
# n1 <- 1000
# n2 <- 1000
# u <- c(rep(p,n1),rep(q,n2))
# lambda <- sum(u^2)
# u <- u/sqrt(lambda)
# B <- matrix(c(p^2, p*q, p*q, q^2), nrow = 2)
# u1 <- B[c(rep(1,n1),rep(2,n2)),c(rep(1,n1),rep(2,n2))] %*% u
# s1 <- c(rep(1,n1),rep(0,n2))
# s2 <- c(rep(0,n1),rep(1,n2))
# s1.tilde <- s1 - u*sum(u*s1)
# s2.tilde <- s2 - u*sum(u*s2)
# # ss1 <- B[c(rep(1,n1),rep(2,n2)),c(rep(1,n1),rep(2,n2))] %*% s1.tilde
# # ss2 <- B[c(rep(1,n1),rep(2,n2)),c(rep(1,n1),rep(2,n2))] %*% s2.tilde
# for(i in 1:nmc){
# A <- sbm.game(n, B, c(n1,n2), loops = TRUE)
# # yy <- irlba(A[],1)
# # Xhat <- as.vector(yy$u %*% as.numeric(sqrt(yy$d)))
# # phat[i] <- mean(Xhat[1:n1])
# # qhat[i] <- mean(Xhat[(n1+1):n])
#  mmu <- A[] %*% u
#  mms1 <- A[] %*% s1.tilde
#  mms2 <- A[] %*% s2.tilde
#  dphat[i] <- sum(mms1*(mmu - u1))/(lambda*sqrt(lambda))
#  dqhat[i] <- sum(mms2*(mmu - u1))/(lambda*sqrt(lambda))
# }
#
# biasp2q2 <- function(rho,p,q){
# delta <- rho*p^2 + (1 - rho)*q^2
# a1 <- 1 - rho*p^2/delta
# b1 <- rho*a1*p/delta^2*(rho*p^2*(1 - p^2) + (1 - rho)*p*q*(1 - p*q))
# b2 <- (1 - rho)*(rho*p*q)*q/delta^3*(rho*p*q*(1 - p*q) + (1 - rho)*q^2*(1 - q^2))
# return(b1 - b2)
# }
#
# plot1 <- function(){
# pseq <- seq(from = 0.01,to = 0.99,by = 0.01)
# qseq <- seq(from = 0.01,to = 0.99,by = 0.01)
#
# error.p <- matrix(0, length(pseq),length(qseq))
# error.q <- matrix(0, length(pseq),length(qseq))
# error.p.ind <- matrix(0, length(pseq),length(qseq))
# error.q.ind <- matrix(0, length(pseq),length(qseq))
# error.p.FIM <- error.p.ind
# error.q.FIM <- error.q.ind
#
# rho <- 0.5
# for(i in 1:length(pseq)){
#   for(j in 1:length(qseq)){
#     p <- pseq[i]
#     q <- qseq[j]
# varp.noindependent <- p2q2(rho,p,q)[1]
# varp.independent <- p2q2(rho,p,q)[2]
# varq.noindependent <- p2q2(1 - rho,q,p)[1]
# varq.independent <- p2q2(1 - rho, q, p)[2]
#
# biasp <- biasp2q2(rho,p,q)
# biasq <- biasp2q2(1 - rho,q,p)
#
# error.p[i,j] <- varp.noindependent
# error.q[i,j] <- varq.noindependent
#
# error.p.ind[i,j] <- varp.independent
# error.q.ind[i,j] <- varq.independent
#
# error.p.FIM[i,j] <- FIM(rho,p,q)[1,1]*rho^2
# error.q.FIM[i,j] <- FIM(1 - rho,q,p)[1,1]*(1 - rho)^2
# }
# }
# levelplot(error.p/error.p.ind,col.regions = new.palette(20), main = "Ratio of variance for p")
# levelplot(error.q/error.q.ind,col.regions = new.palette(20),main = "Ratio of variance for q")
# levelplot(error.p/error.p.FIM,col.regions = new.palette(20), main = "Ratio of variance for p FIM")
# levelplot(error.q/error.q.FIM,col.regions = new.palette(20),main = "Ratio of variance for q FIM")
# }

FIM <- function(rho,p,q){
  a11 <- 2*rho^2/(1 - p^2) + rho*(1 - rho)*q/(p*(1 - p*q))
  a12 <- rho*(1 - rho)/(1 - p*q)
  a22 <- 2*(1 - rho)^2/(1 - q^2) + rho*(1 - rho)*p/(q*(1 - p*q))

  A <- matrix(c(a11,a12,a12,a22),nrow = 2)

  return(solve(A))
}

experiment3.2toinfty <- function(nmc=100, nseq=seq(from = 1000, to = 4000, by = 250)){
  start.time <- Sys.time()
  rho <- c(0.6,0.4)
  B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)
  x <- eigen(B)
  xx <- x$vectors %*% diag(sqrt(x$values))

  all <- data.frame(n=numeric(0),error=numeric(0),name=character(0))
  error.max = 0
  for(n in nseq){
    cat('n=',n,'\n')
    for(j in 1:nmc){
      A <- rg.SBM(n, B, rho,condition=TRUE)
      tau <- A$tau
      X <- xx[tau,]

      Xhat <- stfp(A$adjacency,dim = 2, method="svd")

      W <- procrustes(Xhat,X)
      R <- Xhat%*%W - X
      sq.error <- sqrt(max(rowSums(R^2)))
      all = rbind(all,data.frame(n=n,sq.error=sq.error,name="Observed"))
#      all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(Xhat)),name="Predicted"))
      cat('\rMC = ',j,' ',format(Sys.time()-start.time))
    }
#    all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(X)),name="Expected"))
    cat('\n')

    p<- ggplot(all)+aes(factor(n),sq.error) + geom_violin() + xlab("number of vertices") + ylab("error")
    ggsave("L2_Error_Violin.pdf",width=6,height=4)
  }

  all
}

experiment3.2toinfty.lse <- function(nmc=100, nseq=seq(from = 1000, to = 4000, by = 250)){
  start.time <- Sys.time()
  rho <- c(0.6,0.4)
  B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)
  x <- eigen(B)
  xx <- x$vectors %*% diag(sqrt(x$values))

  all <- data.frame(n=numeric(0),error=numeric(0),name=character(0))
  error.max = 0
  for(n in nseq){
    cat('n=',n,'\n')
    for(j in 1:nmc){
      A <- sbm.game(n, B, rho*n)
      tau <- rep(1:2,rho*n)
      X <- xx[tau,]
      deg <- rep(c(rho[1]*B[1,1] + rho[2]*B[1,2], rho[1]*B[2,1] + rho[2]*B[2,2]),rho*n)*n
      Xtilde <- X/sqrt(deg)

      Xbreve <- embed_laplacian_matrix(A, 2, type = "DAD",options=list(maxiter=100000))$X

      W <- procrustes(Xbreve,Xtilde)
      R <- Xbreve%*%W - Xtilde
      sq.error <- sqrt(max(rowSums(R^2)))
      all = rbind(all,data.frame(n=n,sq.error=sq.error,name="Observed"))
#      all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(Xhat)),name="Predicted"))
      cat('\rMC = ',j,' ',format(Sys.time()-start.time))
    }
#    all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(X)),name="Expected"))
    cat('\n')

    p<- ggplot(all)+aes(factor(n),sq.error) + geom_violin() + xlab("number of vertices") + ylab("error")
    ggsave("L2_Error_Violin.pdf",width=6,height=4)
  }

  all
}



#
# ggplot(all, aes(x = n, y = sq.error)) +
#     geom_errorbar(aes(ymin = sq.error - sd, ymax = sq.error + sd), colour = "black", width = .1,
#                   position = pd) +
#     geom_point(position = pd, size = 3, shape = 21, fill="white") +
#     xlab("number of vertices") +
#         ylab("error in 2 to infty norm") + theme_light() + ylim(0,0.8) +
#     theme(axis.text = element_text(size = 15)) +
#     stat_function(fun = function(x) x^(-0.54)*10^(1.31), linetype = "dashed")
# ggsave("2toinfty_ase.pdf")
#
# ggplot(all.lse, aes(x = n, y = sq.error)) +
#     geom_errorbar(aes(ymin = sq.error - sd, ymax = sq.error + sd), colour = "black", width = .1,
#                   position = pd) +
#     geom_point(position = pd, size = 3, shape = 21, fill="white") +
#     xlab("number of vertices") +
#     ylab("error in 2 to infty norm") + theme_light() + ylim(0,0.04) +
#     theme(axis.text = element_text(size = 15)) +
#     stat_function(fun = function(x) x^(-1.06)*10^(1.57), linetype = "dashed")
# ggsave("2toinfty_lse.pdf")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

experiment1.lse <- function(n = 1000){
    #n <- 8000
    rho <- c(0.6,0.4)
    B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)

    A <- sbm.game(n, B, rho*n)
    tau <- rep(1:length(rho), rho*n)
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))
    X <- xx[tau,]
    deg <- c(rho[1]*B[1,1] + rho[2]*B[1,2], rho[1]*B[2,1] + rho[2]*B[2,2])*n
    deg.vec <- rep(deg, rho*n)
    Xtilde <- X/sqrt(deg)

    Xhat <- embed_laplacian_matrix(A, no = 2, type = "DAD", option = list(maxiter = 10000))$X
    W <- procrustes(Xhat,Xtilde)
    residue <- sqrt(n)*(Xhat%*%W - Xtilde)
    print(norm(Xhat%*%W - Xtilde,"F")^2)

    sigs <- block.var(xx,rho)

    print(sigs)

   #  mean.residue1 <- colMeans(residue[which(tau == 1),])
   sd.residue1 <- sd.multivariate(residue[which(tau == 1),])
   #  mean.residue2 <- colMeans(residue[which(tau == 2),])
   sd.residue2 <- sd.multivariate(residue[which(tau == 2),])
   print(sd.residue1)
   print(sd.residue2)

    df <- data.frame(x = sqrt(n)*(Xhat%*%W)[,1], y = sqrt(n)*(Xhat%*%W)[,2],
                     block = as.character(tau),
                     tau = tau)

    ## library(ellipse)
    ## df_ell <- data.frame()

    ## for(g in unique(df$tau)){
    ##     sig <- sigs[[g]]
    ##     df_ell <- rbind(df_ell, cbind(as.data.frame(
    ##       with(df[df$tau==g,],
    ##            ellipse(sig[1,2]/sqrt(sig[1,1]*sig[2,2]),
    ##            scale=c(sqrt(sig[1,1])/sqrt(n),sqrt(sig[2,2])/sqrt(n)),
    ##            centre=c(xx[g,1],xx[g,2]))),group=g)))

    ## }

    library(ggplot2)
    ggplot(data = df, aes(x = x, y=y,color=block)) + geom_point(size = 2, alpha = 0.3) +
        stat_ellipse(type = "norm", size = 2, linetype = "dashed") +
        geom_point(x = sqrt(n)*xx[1,1]/sqrt(deg[1]),y = sqrt(n)*xx[1,2]/sqrt(deg[1]), size = 6,
                   shape = 16, col = "black") +
        geom_point(x = sqrt(n)*xx[2,1]/sqrt(deg[2]), y = sqrt(n)*xx[2,2]/sqrt(deg[2]), size = 6,
                   shape = 17, col = "black") +
        theme_light() + theme(legend.position = "none")  +
        theme(axis.text = element_text(size = 15))
   ggsave(paste("clusplot",n,"_lse.pdf",sep=''),width=7,height=7)

}

experiment1.ase <- function(n = 1000){
    #n <- 8000
    rho <- c(0.6,0.4)
    B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)

    A <- sbm.game(n, B, rho*n)
    tau <- rep(1:length(rho), rho*n)
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))
    X <- xx[tau,]

    Xhat <- embed_adjacency_matrix(A, no = 2)$X
    W <- procrustes(Xhat,X)
    residue <- sqrt(n)*(Xhat%*%W - X)
    print(norm(Xhat%*%W - X,"F")^2)

    sigs <- block.var(xx,rho)

    print(sigs)

   #  mean.residue1 <- colMeans(residue[which(tau == 1),])
   sd.residue1 <- sd.multivariate(residue[which(tau == 1),])
   #  mean.residue2 <- colMeans(residue[which(tau == 2),])
   sd.residue2 <- sd.multivariate(residue[which(tau == 2),])
   print(sd.residue1)
   print(sd.residue2)

    df <- data.frame(x = (Xhat%*%W)[,1], y = (Xhat%*%W)[,2],
                     block = as.character(tau),
                     tau = tau)

    ## library(ellipse)
    ## df_ell <- data.frame()

    ## for(g in unique(df$tau)){
    ##     sig <- sigs[[g]]
    ##     df_ell <- rbind(df_ell, cbind(as.data.frame(
    ##       with(df[df$tau==g,],
    ##            ellipse(sig[1,2]/sqrt(sig[1,1]*sig[2,2]),
    ##            scale=c(sqrt(sig[1,1])/sqrt(n),sqrt(sig[2,2])/sqrt(n)),
    ##            centre=c(xx[g,1],xx[g,2]))),group=g)))

    ## }

    library(ggplot2)
    ggplot(data = df, aes(x = x, y=y, color=block)) +
      geom_point(size = 2, alpha = 0.3) +
        stat_ellipse(type = "norm", size = 2, linetype = "dashed") #+
        geom_point(x = xx[1,1],y = xx[1,2], size = 6,
                   shape = 16, col = "black") +
        geom_point(x = xx[2,1], y = xx[2,2], size = 6,
                   shape = 17, col = "black") +
        theme_light() + theme(legend.position = "none")  +
        theme(axis.text = element_text(size = 15))
   ggsave(paste("clusplot",n,"_ase.pdf",sep=''),width=7,height=7)

}

