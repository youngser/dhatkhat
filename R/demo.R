## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------------
#knitr::opts_chunk$set(fig.align = 'center', message = FALSE)


## ----pkg-------------------------------------------------------------------------------------------------------------------
## required packages and sources
options(repos=structure(c(CRAN="http://cloud.r-project.org/")))
requiredPkg <- c("tidyverse", "igraph", "Matrix", "mclust")
newPkg <- requiredPkg[!(requiredPkg %in% installed.packages()[, "Package"])]
if (length(newPkg)) install.packages(newPkg, dependencies = TRUE)
update.packages(newPkg, ask = FALSE)
invisible(sapply(requiredPkg, require, character.only = TRUE))

source(url("http://www.cis.jhu.edu/~parky/dhatKhat/R/utils.R"))


## ----param-----------------------------------------------------------------------------------------------------------------
## parameters
n <- 100 # number of vertices
p <- 0.09 # a probability that between blocks can have edges
B <- cbind( c(.2, p), c(p, .1) ) # blockconnectivityprobabilitymatrix
rho <- c(.5, .5) # block membership probabilities, \pi in the paper
K <- length(rho) # number of blocks
Y <- rep(1:K, times=rho*n) # block membership vector, vertex labels


## ----g---------------------------------------------------------------------------------------------------------------------
## generate a graph
set.seed(42+49182)
g <- sbm.game(n, pref.matrix = as.matrix(B),
             block.sizes = rho*n, directed = F, loops = F); summary(g)
#image(as_adjacency_matrix(g), main="Adjacency Matrix")
plotA(g)


## ----emb-------------------------------------------------------------------------------------------------------------------
## embed graph
dmax <- 8 # maximum embedding dimension
ase <- embed.graph(A=g[,,sparse=F], dim=dmax, scaling = TRUE)
elb <- getElbows(abs(ase$D), n=3, plot=TRUE); title("Scree plot with 3 suggested elbows")


## ----mc--------------------------------------------------------------------------------------------------------------------
## run our algorithms
Kmax <- 4 # maximum number of clusters
bic <- ari <- rep(-Inf, dmax-1)
Khat <- rep(0, dmax-1)
out.SMS <- NULL; out.SMS$bic <- -Inf
for (d in 1 : (dmax - 1)) {
  for (k in 1 : Kmax) {
    # SMS_Reduced
    mout.SMS_R <- bicGmmVVV(ase$X[, 1:d, drop=FALSE], k)
    if (mout.SMS_R$bic > bic[d]) {
      bic[d] <- mout.SMS_R$bic
      ari[d] <- adjustedRandIndex(Y, mout.SMS_R$classification)
      Khat[d] <- k
    }

    # SMS
    mout.SMS <- bicGmmBlock2(ase$X, d, k)
    if (mout.SMS$bic > out.SMS$bic) {
      out.SMS <- mout.SMS
    }
  }
}

## calculate ARIs
## SMS
ari_SMS <- adjustedRandIndex(Y, out.SMS$classification)
dhat_SMS <- out.SMS$d
khat_SMS <- out.SMS$k

## SMS_Reduced
ari_SMS_R <- ari[dhat_SMS]
dhat_SMS_R <- dhat_SMS
khat_SMS_R <- khat_SMS

## ZG1
dhat_ZG1 <- min(elb[1], dmax-1)
ari_ZG1 <- ari[dhat_ZG1]
khat_ZG1 <- Khat[dhat_ZG1]

## ZG2
dhat_ZG2 <- min(elb[2], dmax-1)
ari_ZG2 <- ari[dhat_ZG2]
khat_ZG2 <- Khat[dhat_ZG2]

## ZG3
dhat_ZG3 <- min(elb[3], dmax-1)
ari_ZG3 <- ari[dhat_ZG3]
khat_ZG3 <- Khat[dhat_ZG3]


## ----alg-------------------------------------------------------------------------------------------------------------------
## Louvain
result_Louvain <- cluster_louvain(g)
ari_Louvain <- adjustedRandIndex(Y, result_Louvain$membership)
khat_Louvain <- max(membership(result_Louvain))

## Walktrap
result_Rw <- cluster_walktrap(g)
ari_Walktrap <- adjustedRandIndex(Y, result_Rw$membership)
khat_Walktrap <- max(membership(result_Rw))

## IRM
if (!file.exists("RData/irm_out.RData")) {
  system.time(result_IRM <- irm(g[,,sparse=FALSE], sweeps = 1000)) # ~10 min to run on a macbook
  ari_IRM <- max(sapply(1:nrow(result_IRM), function(x) adjustedRandIndex(Y, result_IRM[x,])))
  khat_IRM <- max(mode.irm(result_IRM))
#ari_IRM <- khat_IRM <- NA
  dir.create("RData", showWarnings = FALSE)
  save(result_IRM, ari_IRM, khat_IRM, file="RData/irm_out.RData")
} else {
  load("RData/irm_out.RData")
}


## ----output----------------------------------------------------------------------------------------------------------------
algs <- c("SMS", "SMS_R", "ZG1", "ZG2", "ZG3", "Louvain", "Walktrap", "IRM")
df.out <- data.frame(method = factor(algs, levels=algs),
                     ARI = c(ari_SMS, ari_SMS_R, ari_ZG1, ari_ZG2, ari_ZG3,
                             ari_Louvain, ari_Walktrap, ari_IRM),
                     dhat = c(dhat_SMS, dhat_SMS_R, dhat_ZG1, dhat_ZG2, dhat_ZG3,
                              NA, NA, NA),
                     Khat = c(khat_SMS, khat_SMS_R, khat_ZG1, khat_ZG2, khat_ZG3,
                              khat_Louvain, khat_Walktrap, khat_IRM))
df.out

df.out %>% ggplot(aes(method, ARI, color=method, fill=method)) +
  geom_col(alpha=0.5, show.legend = FALSE) +
  geom_text(aes(label=paste0("hat(K) == ", Khat)), parse=TRUE, vjust=-0.2, show.legend = FALSE) +
  theme(axis.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=12)) +
  #  theme(axis.ticks = element_line(size = 5)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(strip.text=element_text(size=rel(1.2))) +
  theme(legend.text = element_text(colour="black", size = 15, face = "plain"))

