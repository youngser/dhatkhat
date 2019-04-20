# arguments get from command line when running the codes in MARCC
args <-commandArgs(trailingOnly = TRUE)
fileIndex <- as.numeric(args[1])

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
############################### END REQURIES ########################


############################### PARAMETERS ###########################
dmax <- 100
dd <- 50
kmax <- 30
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
ase <- stfp(A = Adj, dim=dmax, scaling = TRUE)
X <- ase$X

labels <- rep(0, dim(Adj)[1])
labels0 <- rep(0, dim(Adj)[1])
labels1 <- rep(0, dim(Adj)[1])
labels2 <- rep(0, dim(Adj)[1])
labels3 <- rep(0, dim(Adj)[1])
###############################END  PARAMETERS ###########################
print(paste("fileIndex=", fileIndex))
for (d in 1 : dd) {
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
khat0 <- mcOut0$G
mcOut1 <- Mclust(ase$X[,1:elb[1]], G = 1 : kmax)
khat1 <- mcOut1$G
mcOut2 <- Mclust(ase$X[,1:elb[2]], G = 1 : kmax)
khat2 <- mcOut2$G
mcOut3 <- Mclust(ase$X[,1:elb[3]], G = 1 : kmax)
khat3 <- mcOut3$G

filename = paste("RS_DS01216_", dataList[fileIndex, 1], ".RData", sep="")

save.image(filename)