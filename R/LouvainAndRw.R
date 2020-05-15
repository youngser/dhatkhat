suppressMessages(require(igraph))
suppressMessages(require(mclust))

source("utils.R")

args <-commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
mc <- as.numeric(args[2])

pList <- c(0.09, 0.095, 0.1, 0.105, 0.11, 0.115)
n <- 500
nmc <- 100
sweeps <- 100
meanAriLouvain <- rep(0, length(pList))
meanAriWalktrap <- rep(0, length(pList))
meanAriIRM <- rep(0, length(pList))

rho <- c(.5, .5)
m <- length(rho)
Y <- rep(1:m, times=rho*n)

for (mc in 1 : nmc) {
  set.seed(mc * 7 + 49182)

  for (i in 1 : length(pList)) {
    p <-pList[i]
    # 2-block SBM
    B <- cbind( c(.2, p), c(p, .1) )

    # # 3-block SBM
    # B <- cbind( c(.2, p, p), c(p, .2, p), c(p, p, .2))
    # rho <- c(.33, .33, .34)
    # m <- length(rho)
    # Y <- rep(1:m, times=rho*n)

    g <- sbm.game(n, pref.matrix = as.matrix(B),
                 block.sizes = rho*n, directed = F, loops = F)
    A <- as.matrix(g[])

    result_louvain <- cluster_louvain(g)
    result_Rw <- cluster_walktrap(g)
    result_irm <- irm(A, sweeps = 1000)

    meanAriLouvain[i] <- meanAriLouvain[i] +
      adjustedRandIndex(Y, result_louvain$membership) / nmc
    meanAriWalktrap[i] <- meanAriWalktrap[i] +
      adjustedRandIndex(Y, result_Rw$membership) / nmc
    meanAriIRM[i] <-meanAriIRM[i] +
      max(sapply(1:nrow(result_irm), function(x) adjustedRandIndex(Y, result_irm[x,]))) / nmc
  }
}

