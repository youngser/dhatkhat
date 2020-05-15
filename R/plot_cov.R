############################### REQURIES ###########################
library("lattice")
source("toolFun.R")
suppressMessages(require(igraph))
suppressMessages(require(Matrix))
############################### END REQURIES ########################
#set.seed(92620)

B <- cbind(c(0.2, 0.1), c(0.1, 0.25))
rho <- c(0.5, 0.5)
m <- length(rho)
dmax <- 20
#pdf('Figure_4.pdf')

n <- 200
#n <- 2000
nmc <- 1
Y <- rep(1:m, times=rho*n)
y <- matrix(0, dmax, dmax)
for (j in 1:nmc) {
  g<- sbm.game(n, pref.matrix = as.matrix(B),
               block.sizes = rho*n, directed = F, loops = F) 
  ase <- eigASE(A=g[,,sparse=F], dim=dmax, scaling = TRUE)
  
  y <- y + cov(ase$X[1:(n/2),]) / nmc
}

levelplot(y)
#dev.off()
