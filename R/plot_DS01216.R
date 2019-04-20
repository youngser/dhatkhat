############################### REQURIES ###########################
source("SMS_EM.R")

requiredPkg <- c("igraph", "Matrix", "mclust", "dplyr", "clustvarsel", "MASS", "lattice", "ggplot2", "googledrive")
newPkg <- requiredPkg[!(requiredPkg %in% installed.packages()[, "Package"])]
if (length(newPkg)) install.packages(newPkg, dependencies = TRUE)
sapply(requiredPkg, require, character.only = TRUE)

############################### END REQURIES ########################
dList <- read.table("./filename.txt")
ngraphs <- dim(dList)[1]
crange <- 12

dhat_MCG <- dhat_ZG1 <- dhat_ZG2 <- dhat_ZG3 <- rep(0, ngraphs)
khat_MCG <- khat_MCEG <- khat_ZG1 <- khat_ZG2 <- khat_ZG3 <- rep(0, ngraphs)
ari_MCG_tissue <- ari_MCG_Y <- ari_MCG_region <- ari_MCG_hemisphere <- rep(0, ngraphs)
ari_MCEG_tissue <- ari_MCEG_Y <- ari_MCEG_region <- ari_MCEG_hemisphere <- rep(0, ngraphs)
ari_ZG1_tissue <- ari_ZG1_Y <- ari_ZG1_region <- ari_ZG1_hemisphere <- rep(0, ngraphs)
ari_ZG2_tissue <- ari_ZG2_Y <- ari_ZG2_region <- ari_ZG2_hemisphere <- rep(0, ngraphs)
ari_ZG3_tissue <- ari_ZG3_Y <- ari_ZG3_region <- ari_ZG3_hemisphere <- rep(0, ngraphs)
ari_Louvain_tissue <- ari_Louvain_Y <- ari_Louvain_region <- ari_Louvain_hemisphere <- rep(0, ngraphs)
ari_Rw_tissue <- ari_Rw_Y <- ari_Rw_region <- ari_Rw_hemisphere <- rep(0, ngraphs)

dmax <- 100
maxBic <- rep(-1e10, ngraphs)
dd <- 50
kmax <- 30

drive_auth_config(active = FALSE)
folder_url <- "https://drive.google.com/open?id=1zBuIMdpNQFoWZvh7JOepswwFCq2l9Wco"
folder <- drive_get(as_id(folder_url))
files <- drive_ls(folder)


for (fIndex in 1 : ngraphs) {
  print(fIndex)
  filename <- paste("RS_DS01216_", dList[fIndex, 1], ".RData", sep = "")
  filei <- files %>% filter(name == filename)
  file <- drive_download(filei, overwrite = TRUE, verbose = FALSE)
  load(file$local_path)

  bic <- matrix(rep(0, dd*kmax), dd,kmax)
  for (d in 1 : dd) {
    for (k in 1 : kmax) {
      bic[d, k] <- moutBlock2[[d]][[k]]$bic
    }
  }
  c = 1
  nn <- dim(X)[1]
  dkPairs <- maxBIC1(bic, globalMax = 0, regression = 2, penalty = c, nn, D)
  dhat_MCG[fIndex] <- dkPairs$dhat
  khat_MCG[fIndex] <- dkPairs$khat
  
  result_louvain <- cluster_louvain(g2)
  result_Rw <- cluster_walktrap(g2)
  ari_Louvain_tissue[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                                 result_louvain$membership)
  ari_Louvain_hemisphere[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                                     result_louvain$membership)
  ari_Louvain_Y[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                            result_louvain$membership)
  ari_Louvain_region[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                                 result_louvain$membership)
  ari_Rw_tissue[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                            result_Rw$membership)
  ari_Rw_hemisphere[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                                result_Rw$membership)
  ari_Rw_Y[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                       result_Rw$membership)
  ari_Rw_region[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                            result_Rw$membership)
  ari_MCG_tissue[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                      moutBlock2[[dhat_MCG[fIndex]]][[khat_MCG[fIndex]]]$classification)
  ari_MCG_Y[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                 moutBlock2[[dhat_MCG[fIndex]]][[khat_MCG[fIndex]]]$classification)
  ari_MCG_region[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                      moutBlock2[[dhat_MCG[fIndex]]][[khat_MCG[fIndex]]]$classification)
  ari_MCG_hemisphere[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                          moutBlock2[[dhat_MCG[fIndex]]][[khat_MCG[fIndex]]]$classification)
  ari_MCEG_tissue[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                       mcOut0$classification)
  ari_MCEG_Y[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                  mcOut0$classification)
  ari_MCEG_region[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                       mcOut0$classification)
  ari_MCEG_hemisphere[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                           mcOut0$classification)
  ari_ZG1_tissue[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                      mcOut1$classification)
  ari_ZG1_Y[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                 mcOut1$classification)
  ari_ZG1_region[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                      mcOut1$classification)
  ari_ZG1_hemisphere[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                          mcOut1$classification)
  ari_ZG2_tissue[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                      mcOut2$classification)
  ari_ZG2_Y[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                 mcOut2$classification)
  ari_ZG2_region[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                      mcOut2$classification)
  ari_ZG2_hemisphere[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                          mcOut2$classification)
  ari_ZG3_tissue[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "tissue"),
                                      mcOut3$classification)
  ari_ZG3_Y[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "Y"),
                                 mcOut3$classification)
  ari_ZG3_region[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "region"),
                                      mcOut3$classification)
  ari_ZG3_hemisphere[fIndex] <- adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
                                          mcOut3$classification)
  dhat_ZG1[fIndex] <- elb[1]
  dhat_ZG2[fIndex] <- elb[2]
  dhat_ZG3[fIndex] <- elb[3]
  khat_ZG1[fIndex] <- khat1
  khat_ZG2[fIndex] <- khat2
  khat_ZG3[fIndex] <- khat3
  khat_MCEG <mcOut0$G
  
  system(paste("rm ", file$local_path))
}

# plot figure 8
df <- data.frame(ari.sms=ari_MCG_Y, ari_ZG1_Y, ari_ZG2_Y, ari_ZG3_Y,
                 dhat_MCG, dhat_ZG1, dhat_ZG2, dhat_ZG3,
                 khat_MCG, khat_ZG1, khat_ZG2, khat_ZG3)
df2 <- rbind(data.frame(dhat=dhat_MCG, Khat=khat_MCG, ari=ari_MCG_Y, method="SMS"),
             data.frame(dhat=dhat_ZG1, Khat=khat_ZG1, ari=ari_ZG1_Y, method="ZG1"),
             data.frame(dhat=dhat_ZG2, Khat=khat_ZG2, ari=ari_ZG2_Y, method="ZG2"),
             data.frame(dhat=dhat_ZG3, Khat=khat_ZG3, ari=ari_ZG3_Y, method="ZG3"))

suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(latex2exp))

set.seed(123)
library(RColorBrewer)
mycol <- rev(colorRampPalette(brewer.pal(9,"RdYlBu"))(11))
p2 <- ggplot(df2, aes(x=dhat, y=Khat, fill=as.numeric(as.character(ari)))) +
  geom_jitter(aes(shape=method, size=ari), width=0.5, height=0.5, alpha=1) + 
  scale_fill_gradientn(colours=mycol) +
  #       scale_fill_gradient2(low=mycol[1], mid=mycol[5], high=mycol[11],
  # 					 midpoint=0.2,
  # 					 breaks=c(0.01, 0.20, 0.4) ) +
  scale_shape_manual(values=c(21:24)) +
  scale_size(guide = 'none') +
  labs(x=TeX('$\\hat{d}$'), y=TeX('$\\hat{K}$'), color = "ARI", fill="ARI")
p2

# plot Figure 9
hist(ari_MCG_hemisphere - ari_ZG1_hemisphere, xlab = "ARI_MCG - ARI_BIC.ZG1", 
     main="", col = "blue", density = 50)

hist(ari_MCG_hemisphere - ari_ZG2_hemisphere, xlab = "ARI_MCG - ARI_BIC.ZG2", 
     main="", col = "black", density = 50)

hist(ari_MCG_hemisphere - ari_ZG3_hemisphere, xlab = "ARI_MCG - ARI_BIC.ZG3", 
     main="", col = "purple", density = 50)

hist(ari_MCEG_hemisphere - ari_ZG1_hemisphere, xlab = "ARI_MCEG - ARI_BIC.ZG1", 
     main="", col = "blue", angle = 135, density = 10)

hist(ari_MCEG_hemisphere - ari_ZG2_hemisphere, xlab = "ARI_MCEG - ARI_BIC.ZG2", 
     main="", col = "black", angle = 135, density = 10)

hist(ari_MCEG_hemisphere - ari_ZG3_hemisphere, xlab = "ARI_MCEG - ARI_BIC.ZG3", 
     main="", col = "purple", angle = 135, density = 10)

# the results of Table 10
win_MCG_ZG1_hemisphere <- sum(ari_MCG_hemisphere > ari_ZG1_hemisphere)
pvalue_MCG_ZG1_hemisphere <- 1 - pbinom(win_MCG_ZG1_hemisphere, ngraphs, 0.5)
win_MCG_ZG2_hemisphere <- sum(ari_MCG_hemisphere > ari_ZG2_hemisphere)
pvalue_MCG_ZG2_hemisphere <- 1 - pbinom(win_MCG_ZG2_hemisphere, ngraphs, 0.5)
win_MCG_ZG3_hemisphere <- sum(ari_MCG_hemisphere > ari_ZG3_hemisphere)
pvalue_MCG_ZG3_hemisphere <- 1 - pbinom(win_MCG_ZG3_hemisphere, ngraphs, 0.5)
win_MCG_Louvain_hemisphere <- sum(ari_MCG_hemisphere > ari_Louvain_hemisphere)
pvalue_MCG_Louvain_hemisphere <- 1 - pbinom(win_MCG_Louvain_hemisphere, ngraphs, 0.5)
win_MCG_Rw_hemisphere <- sum(ari_MCG_hemisphere > ari_Rw_hemisphere)
pvalue_MCG_Rw_hemisphere <- 1 - pbinom(win_MCG_Rw_hemisphere, ngraphs, 0.5)

win_MCEG_ZG1_hemisphere <- sum(ari_MCEG_hemisphere > ari_ZG1_hemisphere)
pvalue_MCEG_ZG1_hemisphere <- 1 - pbinom(win_MCEG_ZG1_hemisphere, ngraphs, 0.5)
win_MCEG_ZG2_hemisphere <- sum(ari_MCEG_hemisphere > ari_ZG2_hemisphere)
pvalue_MCEG_ZG2_hemisphere <- 1 - pbinom(win_MCEG_ZG2_hemisphere, ngraphs, 0.5)
win_MCEG_ZG3_hemisphere <- sum(ari_MCEG_hemisphere > ari_ZG3_hemisphere)
pvalue_MCEG_ZG3_hemisphere <- 1 - pbinom(win_MCEG_ZG3_hemisphere, ngraphs, 0.5)
win_MCEG_Louvain_hemisphere <- sum(ari_MCEG_hemisphere > ari_Louvain_hemisphere)
pvalue_MCEG_Louvain_hemisphere <- 1 - pbinom(win_MCEG_Louvain_hemisphere, ngraphs, 0.5)
win_MCEG_Rw_hemisphere <- sum(ari_MCEG_hemisphere > ari_Rw_hemisphere)
pvalue_MCEG_Rw_hemisphere <- 1 - pbinom(win_MCEG_Rw_hemisphere, ngraphs, 0.5)

win_MCG_ZG1_tissue <- sum(ari_MCG_tissue > ari_ZG1_tissue)
pvalue_MCG_ZG1_tissue <- 1 - pbinom(win_MCG_ZG1_tissue, ngraphs, 0.5)
win_MCG_ZG2_tissue <- sum(ari_MCG_tissue > ari_ZG2_tissue)
pvalue_MCG_ZG2_tissue <- 1 - pbinom(win_MCG_ZG2_tissue, ngraphs, 0.5)
win_MCG_ZG3_tissue <- sum(ari_MCG_tissue > ari_ZG3_tissue)
pvalue_MCG_ZG3_tissue <- 1 - pbinom(win_MCG_ZG3_tissue, ngraphs, 0.5)
win_MCG_Louvain_tissue <- sum(ari_MCG_tissue > ari_Louvain_tissue)
pvalue_MCG_Louvain_tissue <- 1 - pbinom(win_MCG_Louvain_tissue, ngraphs, 0.5)
win_MCG_Rw_tissue <- sum(ari_MCG_tissue > ari_Rw_tissue)
pvalue_MCG_Rw_tissue <- 1 - pbinom(win_MCG_Rw_tissue, ngraphs, 0.5)

win_MCEG_ZG1_tissue <- sum(ari_MCEG_tissue > ari_ZG1_tissue)
pvalue_MCEG_ZG1_tissue <- 1 - pbinom(win_MCEG_ZG1_tissue, ngraphs, 0.5)
win_MCEG_ZG2_tissue <- sum(ari_MCEG_tissue > ari_ZG2_tissue)
pvalue_MCEG_ZG2_tissue <- 1 - pbinom(win_MCEG_ZG2_tissue, ngraphs, 0.5)
win_MCEG_ZG3_tissue <- sum(ari_MCEG_tissue > ari_ZG3_tissue)
pvalue_MCEG_ZG3_tissue <- 1 - pbinom(win_MCEG_ZG3_tissue, ngraphs, 0.5)
win_MCEG_Louvain_tissue <- sum(ari_MCEG_tissue > ari_Louvain_tissue)
pvalue_MCEG_Louvain_tissue <- 1 - pbinom(win_MCEG_Louvain_tissue, ngraphs, 0.5)
win_MCEG_Rw_tissue <- sum(ari_MCEG_tissue > ari_Rw_tissue)
pvalue_MCEG_Rw_tissue <- 1 - pbinom(win_MCEG_Rw_tissue, ngraphs, 0.5)

win_MCG_ZG1_Y <- sum(ari_MCG_Y > ari_ZG1_Y)
pvalue_MCG_ZG1_Y <- 1 - pbinom(win_MCG_ZG1_Y, ngraphs, 0.5)
win_MCG_ZG2_Y <- sum(ari_MCG_Y > ari_ZG2_Y)
pvalue_MCG_ZG2_Y <- 1 - pbinom(win_MCG_ZG2_Y, ngraphs, 0.5)
win_MCG_ZG3_Y <- sum(ari_MCG_Y > ari_ZG3_Y)
pvalue_MCG_ZG3_Y <- 1 - pbinom(win_MCG_ZG3_Y, ngraphs, 0.5)
win_MCG_Louvain_Y <- sum(ari_MCG_Y > ari_Louvain_Y)
pvalue_MCG_Louvain_Y <- 1 - pbinom(win_MCG_Louvain_Y, ngraphs, 0.5)
win_MCG_Rw_Y <- sum(ari_MCG_Y > ari_Rw_Y)
pvalue_MCG_Rw_Y <- 1 - pbinom(win_MCG_Rw_Y, ngraphs, 0.5)

win_MCEG_ZG1_Y <- sum(ari_MCEG_Y > ari_ZG1_Y)
pvalue_MCEG_ZG1_Y <- 1 - pbinom(win_MCEG_ZG1_Y, ngraphs, 0.5)
win_MCEG_ZG2_Y <- sum(ari_MCEG_Y > ari_ZG2_Y)
pvalue_MCEG_ZG2_Y <- 1 - pbinom(win_MCEG_ZG2_Y, ngraphs, 0.5)
win_MCEG_ZG3_Y <- sum(ari_MCEG_Y > ari_ZG3_Y)
pvalue_MCEG_ZG3_Y <- 1 - pbinom(win_MCEG_ZG3_Y, ngraphs, 0.5)
win_MCEG_Louvain_Y <- sum(ari_MCEG_Y > ari_Louvain_Y)
pvalue_MCEG_Louvain_Y <- 1 - pbinom(win_MCEG_Louvain_Y, ngraphs, 0.5)
win_MCEG_Rw_Y <- sum(ari_MCEG_Y > ari_Rw_Y)
pvalue_MCEG_Rw_Y <- 1 - pbinom(win_MCEG_Rw_Y, ngraphs, 0.5)
