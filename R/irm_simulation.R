# arguments get from command line when running the codes in MARCC
args <-commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
mc <- as.numeric(args[2])

############################### REQURIES ###########################
source("ModelSelect.R")
source("SMS_EM.R")
source("toolFun.R")
source("utils.R")
suppressMessages(require(igraph))
suppressMessages(require(Matrix))
suppressMessages(require(mclust))
suppressMessages(require(dplyr))
#suppressMessages(require(clustvarsel))
############################### END REQURIES ########################

############################### PARAMETERS ###########################
pList <- c(0.09, 0.095, 0.1, 0.105, 0.11, 0.115)
pNameList <- c("009", "0095", "01", "0105", "011", "0115")
p <- pList[i]
set.seed(mc * 7 + 49182)
n <- 500
BF1 <- cbind( c(.2, p), c(p, .1) )
#BF2 <- cbind( c(.2, p, p), c(p, .2, p), c(p, p, .2))

# 2-block SBM
B <- BF1
rho <- c(.5, .5)

## 3-block SBM
# B <- BF2
# rho <- c(.33, .33, .34)

m <- length(rho)
Y <- rep(1:m, times=rho*n)
dmax <- 8
Kmax <- 4
###############################END  PARAMETERS ###########################
g<- sbm.game(n, pref.matrix = as.matrix(B),
             block.sizes = rho*n, directed = F, loops = F)
#ase <- stfp(A=g[,,sparse=F], dim=dmax, scaling = TRUE)
#ase <- eigASE(A=g[,,sparse=F], dim=dmax, scaling = TRUE)
A <- as.matrix(g[])
elapsed <- system.time(Z <- irm(A, sweeps = 1000))[3]
# user   system  elapsed
# 2972.262    0.812 2973.196

save(Z, file=paste0("Z-i","-mc",mc,".Rbin"))
#top.n(Z)
(ari.mat <- sapply(1:nrow(Z), function(x) adjustedRandIndex(Y, Z[x,])))

save.image(paste0("Zimage-i",i,"-mc",mc,".RData"))


# bic <- rep(-1e10, dmax - 1)
# ari <- rep(-1e10, dmax - 1)
# khat <- rep(0, dmax - 1)
# bicDiag1 <- rep(-1e10, dmax - 1)
# bicDiag2 <- rep(-1e10, dmax - 1)
# bicDiag3 <- rep(-1e10, dmax - 1)
# bicDiag4 <- rep(-1e10, dmax - 1)
# outBlock1 <- list(bic = -1e10, classification = rep(0, n), d = 0, k = 0)
# outBlock2 <- list(bic = -1e10, classification = rep(0, n), d = 0, k = 0)
# outBlock3 <- list(bic = -1e10, classification = rep(0, n), d = 0, k = 0)
# outBlock4 <- list(bic = -1e10, classification = rep(0, n), d = 0, k = 0)
#
# for (d in 1 : (dmax - 1)) {
#   print(paste("d=", d))
#   for (k in 1 : Kmax) {
#     print(paste("k=", k))
#     mout <- bicGmmVVV(ase$X[, 1 : d, drop = FALSE], k)
#     if (mout$bic > bic[d]) {
#       bic[d] <- mout$bic
#       ari[d] <- adjustedRandIndex(Y, mout$classification)
#       khat[d] <- k
#     }
#
#     mout <- bicGmmBlock2(ase$X, d, k)
#     if (mout$bic > outBlock2$bic) {
#       outBlock2 <- mout
#     }
#   }
# }
#
# dhat_Opt <- which.max(ari)
# ari_Opt <- max(ari)
# khat_Opt <- khat[dhat_Opt]
#
# ari_MCG <- adjustedRandIndex(Y, outBlock2$classification)
# dhat_MCG <- outBlock2$d
# ari_MCEG <- ari[dhat_MCG]
# khat_MCG <- outBlock2$k
#
# elb <- getElbows(abs(ase$D),n=3,plot=FALSE)
# dhat_ZG1 <- min(elb[1], dmax - 1)
# ari_ZG1 <- ari[dhat_ZG1]
# khat_ZG1 <- khat[dhat_ZG1]
#
# dhat_ZG2 <- min(elb[2], dmax - 1)
# ari_ZG2 <- ari[dhat_ZG2]
# khat_ZG2 <- khat[dhat_ZG2]
#
# dhat_ZG3 <- min(elb[3], dmax - 1)
# ari_ZG3 <- ari[dhat_ZG3]
# khat_ZG3 <- khat[dhat_ZG3]
#
# filename = paste("eig_BF1_n500_p", pNameList[i], "_mc", mc, ".RData", sep="")
#
# save.image(filename)
