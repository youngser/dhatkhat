fileIndex <- 1
D <- 100
dmax <- 3 # 50
kmax <- 2 # 30

############################### REQURIES ###########################
source("ModelSelect.R")
source("SMS_EM.R")
source("toolFun.R")
suppressMessages(require(igraph))
suppressMessages(require(Matrix))
suppressMessages(require(mclust))
suppressMessages(require(dplyr))
suppressMessages(require(clustvarsel))
suppressMessages(require(MASS))
suppressMessages(require(lattice))
suppressMessages(require(mvtnorm))
suppressMessages(require(reshape2))
suppressMessages(require(ggplot2))
############################### END REQURIES ########################


############################### PARAMETERS ###########################
mout <- NULL
moutBlock2 <- NULL
optMout <- NULL
maxBic <- 1e-10

dataList <- read.table("./filename.txt")
index <- paste("http://www.cis.jhu.edu/~parky/Microsoft/JHU-MSR/ZMx2/BNU1/DS01216/"
               , dataList[fileIndex, 1], ".graphml", sep = "")
g <- read_graph(file = index, format = "graphml")
subInx <- ((vertex_attr(g, "tissue")=="white") | (vertex_attr(g, "tissue")=="gray")) & ((vertex_attr(g, "hemisphere")=="left") | (vertex_attr(g, "hemisphere")=="right"))
g1 <- induced_subgraph(g, subInx)
clu <- components(g1)
cInx <- which.max(clu$csize)
g2 <- induced_subgraph(g1, names(clu$membership[clu$membership==cInx]))
Adj <- as_adjacency_matrix(g2, names = FALSE, sparse = FALSE)
ase <- stfp(A = Adj, dim=D, scaling = TRUE)
#ase <- eigASE(A = Adj, dim=D, scaling = TRUE)
X <- ase$X

labels <- rep(0, dim(Adj)[1])
labels0 <- rep(0, dim(Adj)[1])
labels1 <- rep(0, dim(Adj)[1])
labels2 <- rep(0, dim(Adj)[1])
labels3 <- rep(0, dim(Adj)[1])
###############################END  PARAMETERS ###########################
print(paste("fileIndex=", fileIndex))
for (d in 1 : dmax) {
  print(paste("d=", d))
  for (k in 1 : kmax) {
    print(paste("k=", k))
    mout0 <- Mclust(X[,1:d], G = k:k, prior = priorControl())
    mout[[k]] <- bicGmmBlock2(X, d, k, init = mout0$classification)
    if (mout[[k]]$bic > maxBic) {
      maxBic <- mout[[k]]$bic
      dhat <- d
      khat <- k
      optMout <- mout[[k]]
    }
  }
  moutBlock2[[d]] <- mout
}

elb <- getElbows(abs(ase$D),n=3,plot=FALSE)
mcOut0 <- Mclust(ase$X[,1:dhat], G = 1 : kmax)
khat_MCEG <- mcOut0$G
mcOut1 <- Mclust(ase$X[,1:elb[1]], G = 1 : kmax)
khat_ZG1 <- mcOut1$G
mcOut2 <- Mclust(ase$X[,1:elb[2]], G = 1 : kmax)
khat_ZG2 <- mcOut2$G
mcOut3 <- Mclust(ase$X[,1:elb[3]], G = 1 : kmax)
khat_ZG3 <- mcOut3$G
dhat_ZG1 <- elb[1]
dhat_ZG2 <- elb[2]
dhat_ZG3 <- elb[3]

# computing (d,K) pair by regression
bic <- matrix(rep(0, dmax*kmax), dmax,kmax)
for (d in 1 : dmax) {
  for (k in 1 : kmax) {
    bic[d, k] <- moutBlock2[[d]][[k]]$bic
  }
}
c = 1
nn <- dim(X)[1]
dkPairs <- maxBIC1(bic, globalMax = 0, regression = 2, penalty = c, nn, D)
dhat_MCG <- dkPairs$dhat
khat_MCG <- dkPairs$khat

dhat <- c(dhat_ZG1, dhat_ZG2, dhat_ZG3, dhat_MCG)
Khat <- c(khat_ZG1, khat_ZG2, khat_ZG3, dhat_MCG)
df.dk <- data.frame(
	dhat, Khat, type=c("ZG1","ZG2","ZG3","MCG")
)
df.dk


# evaluate the ARI by 4 types of groud truths
ari_MCG_tissue <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                           moutBlock2[[dhat_MCG]][[khat_MCG]]$classification)
ari_MCG_Y <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                      moutBlock2[[dhat_MCG]][[khat_MCG]]$classification)
ari_MCG_region <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                           moutBlock2[[dhat_MCG]][[khat_MCG]]$classification)
ari_MCG_hemisphere <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                               moutBlock2[[dhat_MCG]][[khat_MCG]]$classification)
ari_MCEG_tissue <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                         mcOut0$classification)
ari_MCEG_Y <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                    mcOut0$classification)
ari_MCEG_region <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                         mcOut0$classification)
ari_MCEG_hemisphere <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                             mcOut0$classification)
ari_ZG1_tissue <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                         mcOut1$classification)
ari_ZG1_Y <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                    mcOut1$classification)
ari_ZG1_region <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                         mcOut1$classification)
ari_ZG1_hemisphere <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                             mcOut1$classification)
ari_ZG2_tissue <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                         mcOut2$classification)
ari_ZG2_Y <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                    mcOut2$classification)
ari_ZG2_region <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                         mcOut2$classification)
ari_ZG2_hemisphere <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                             mcOut2$classification)
ari_ZG3_tissue <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                         mcOut3$classification)
ari_ZG3_Y <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                    mcOut3$classification)
ari_ZG3_region <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                         mcOut3$classification)
ari_ZG3_hemisphere <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                             mcOut3$classification)

tissue <- c(ari_ZG1_tissue,ari_ZG2_tissue, ari_ZG3_tissue, ari_MCEG_tissue, ari_MCG_tissue)
hemisphere <- c(ari_ZG1_hemisphere,ari_ZG2_hemisphere, ari_ZG3_hemisphere, ari_MCEG_hemisphere, ari_MCG_hemisphere)
region <- c(ari_ZG1_region,ari_ZG2_region, ari_ZG3_region, ari_MCEG_region, ari_MCG_region)
Y <- c(ari_ZG1_Y,ari_ZG2_Y, ari_ZG3_Y, ari_MCEG_Y, ari_MCG_Y)

df.ari <- data.frame(
	tissue=tissue,
	hemisphere=hemisphere,
	region=region,
	Y=Y,
	type=c("ZG1","ZG2","ZG3","MCEG","MCG")
)
df.ari
df.ari %>% melt %>%
	ggplot(aes(x=variable, y=value, fill=type, color=type)) +
	geom_col(position = position_dodge(), alpha=0.5) +
	labs(x="label",y="ari")
