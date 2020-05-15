############################### REQURIES ###########################
source("SMS_EM.R")
source("irm.R")

requiredPkg <- c("igraph", "Matrix", "mclust", "dplyr", "clustvarsel", "MASS", "lattice", "ggplot2", "rvest")
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
ari_IRM_tissue <- ari_IRM_Y <- ari_IRM_region <- ari_IRM_hemisphere <- rep(0, ngraphs)

dmax <- 100
maxBic <- rep(-1e10, ngraphs)
dd <- 50
kmax <- 30

# drive_auth_config(active = FALSE)
# folder_url <- "https://drive.google.com/open?id=1zBuIMdpNQFoWZvh7JOepswwFCq2l9Wco"
# folder <- drive_get(as_id(folder_url))
# files <- drive_ls(folder)

folder_url <- "http://www.cis.jhu.edu/~parky/dhatKhat/Results/result_01216_dmax100_d50_181212/"
library(rvest) # the new package, version 0.3.0
docs <- read_html(folder_url)
links <- html_table(docs, fill = TRUE)[[1]]$Name
dfiles <- links[grepl(glob2rx("*.RData"), links)]; length(dfiles)
files <- paste0(folder_url, dfiles)

for (fIndex in 1 : ngraphs) {
#  print(fIndex)
  filename <- paste("RS_DS01216_", dList[fIndex, 1], ".RData", sep = "")
  # filei <- files %>% filter(name == filename)
  # file <- drive_download(filei, overwrite = TRUE, verbose = FALSE)
  # load(file$local_path)
  ufile <- url(files[dfiles %in% filename])
#  ufile <- paste0("../Results/result_01216_dmax100_d50_181212/",files[dfiles %in% filename])
  load(ufile);

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
  # result_irm <- irm(as.matrix(g2[]), sweeps = 100)
  #
  # ari_IRM_tissue[fIndex] <- max(sapply(1:nrow(result_irm), function(x) adjustedRandIndex(vertex_attr(g2, name = "tissue"),
  #                                                                                        result_irm[x,])))
  # ari_IRM_hemisphere[fIndex] <- max(sapply(1:nrow(result_irm), function(x) adjustedRandIndex(vertex_attr(g2, name = "hemisphere"),
  #                                                                                        result_irm[x,])))
  # ari_IRM_Y[fIndex] <- max(sapply(1:nrow(result_irm), function(x) adjustedRandIndex(vertex_attr(g2, name = "Y"),
  #                                                                                        result_irm[x,])))
  # ari_IRM_region[fIndex] <- max(sapply(1:nrow(result_irm), function(x) adjustedRandIndex(vertex_attr(g2, name = "region"),
  #                                                                                        result_irm[x,])))
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
#  khat_MCEG <mcOut0$G

  close(ufile)
}
#save.image(file="DS01216.RData")
print(load("RData/DS01216.RData"))

# plot figure 12
df <- data.frame(ari.sms=ari_MCG_Y, ari.sms.r=ari_MCEG_Y,
                 ari_ZG1_Y, ari_ZG2_Y, ari_ZG3_Y,
                 dhat_MCG, dhat_ZG1, dhat_ZG2, dhat_ZG3,
                 khat_MCG, khat_ZG1, khat_ZG2, khat_ZG3)
df2 <- rbind(data.frame(dhat=dhat_MCG, Khat=khat_MCG, ari=ari_MCG_Y, method="SMS"),
             data.frame(dhat=dhat_ZG1, Khat=khat_ZG1, ari=ari_ZG1_Y, method="ZG1"),
             data.frame(dhat=dhat_ZG2, Khat=khat_ZG2, ari=ari_ZG2_Y, method="ZG2"),
             data.frame(dhat=dhat_ZG3, Khat=khat_ZG3, ari=ari_ZG3_Y, method="ZG3"))
#save(df, df2, file="df-drisophila-sms-ZG.RData")
print(load("RData/df-drisophila-sms-ZG.RData"))

suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(latex2exp))

set.seed(123)
library(RColorBrewer)
mycol <- rev(colorRampPalette(brewer.pal(9,"RdYlBu"))(11))
p2 <- ggplot(df2, aes(x=dhat, y=Khat)) +
	geom_jitter(aes(shape=method, color=method, alpha=ari, size=ari), width=0.5,height=0.5) +
	#	facet_wrap(~method, nrow=2) +
	#	geom_point(aes(size=ari, color=method))
	#	scale_color_gradientn(colours=matlab.like(10)) +
	scale_shape_manual(values=c(16, 17, 15, 18))+
    scale_alpha(name = "ARI") +
    scale_size_continuous(name = "ARI") +
    #    scale_shape_manual(values=c(21:24)) +
	#	scale_size(guide = 'none') +
	labs(x=TeX('$\\hat{d}$'), y=TeX('$\\hat{K}$')) +
#    theme(legend.key.size = unit(3,"point")) +
    theme(axis.title=element_text(size=15)) +
    theme(axis.text.x=element_text(size=15)) +
    #  theme(axis.ticks = element_line(size = 5)) +
    theme(axis.text.y=element_text(size=15)) +
#    theme(strip.text=element_text(size=rel(1.2))) +
    theme(legend.text = element_text(colour="black", size = 15, face = "plain"))
p2 + guides(shape = guide_legend(override.aes = list(size=3)))
ggsave("connectome_dhat_khat.pdf", width=6, height=5)

# plot Figure 13
pdf("connectome_ARI_mcg_zg1.pdf")
hist(ari_MCG_hemisphere - ari_ZG1_hemisphere, xlab = "ARI_SMS - ARI_BIC.ZG1",
     main="", col = "blue", density = 30);
abline(v=0, col=2, lty=2, lwd=5);
dev.off()

pdf("connectome_ARI_mcg_zg2.pdf")
hist(ari_MCG_hemisphere - ari_ZG2_hemisphere, xlab = "ARI_SMS - ARI_BIC.ZG2",
     main="", col = "black", density = 30);
abline(v=0, col=2, lty=2, lwd=5);
dev.off()

pdf("connectome_ARI_mcg_zg3.pdf")
hist(ari_MCG_hemisphere - ari_ZG3_hemisphere, xlab = "ARI_SMS - ARI_BIC.ZG3",
     main="", col = "purple", density = 30);
abline(v=0, col=2, lty=2, lwd=5);
dev.off()

pdf("connectome_ARI_mceg_zg1.pdf")
hist(ari_MCEG_hemisphere - ari_ZG1_hemisphere, xlab = "ARI_SMS_Reduced - ARI_BIC.ZG1",
     main="", col = "blue", angle = 135, density = 10);
abline(v=0, col=2, lty=2, lwd=5);
dev.off()

pdf("connectome_ARI_mceg_zg2.pdf")
hist(ari_MCEG_hemisphere - ari_ZG2_hemisphere, xlab = "ARI_SMS_Reduced - ARI_BIC.ZG2",
     main="", col = "black", angle = 135, density = 10);
abline(v=0, col=2, lty=2, lwd=5);
dev.off()

pdf("connectome_ARI_mceg_zg3.pdf")
hist(ari_MCEG_hemisphere - ari_ZG3_hemisphere, xlab = "ARI_SMS_Reduced - ARI_BIC.ZG3",
     main="", col = "purple", angle = 135, density = 10);
abline(v=0, col=2, lty=2, lwd=5);
dev.off()

# the results of Table 1
(win_MCG_ZG1_hemisphere <- sum(ari_MCG_hemisphere > ari_ZG1_hemisphere))
(pvalue_MCG_ZG1_hemisphere <- 1 - pbinom(win_MCG_ZG1_hemisphere, ngraphs, 0.5))
(win_MCG_ZG2_hemisphere <- sum(ari_MCG_hemisphere > ari_ZG2_hemisphere))
(pvalue_MCG_ZG2_hemisphere <- 1 - pbinom(win_MCG_ZG2_hemisphere, ngraphs, 0.5))
(win_MCG_ZG3_hemisphere <- sum(ari_MCG_hemisphere > ari_ZG3_hemisphere))
(pvalue_MCG_ZG3_hemisphere <- 1 - pbinom(win_MCG_ZG3_hemisphere, ngraphs, 0.5))
(win_MCG_Louvain_hemisphere <- sum(ari_MCG_hemisphere > ari_Louvain_hemisphere))
(pvalue_MCG_Louvain_hemisphere <- 1 - pbinom(win_MCG_Louvain_hemisphere, ngraphs, 0.5))
(win_MCG_Rw_hemisphere <- sum(ari_MCG_hemisphere > ari_Rw_hemisphere))
(pvalue_MCG_Rw_hemisphere <- 1 - pbinom(win_MCG_Rw_hemisphere, ngraphs, 0.5))

(win_MCEG_ZG1_hemisphere <- sum(ari_MCEG_hemisphere > ari_ZG1_hemisphere))
(pvalue_MCEG_ZG1_hemisphere <- 1 - pbinom(win_MCEG_ZG1_hemisphere, ngraphs, 0.5))
(win_MCEG_ZG2_hemisphere <- sum(ari_MCEG_hemisphere > ari_ZG2_hemisphere))
(pvalue_MCEG_ZG2_hemisphere <- 1 - pbinom(win_MCEG_ZG2_hemisphere, ngraphs, 0.5))
(win_MCEG_ZG3_hemisphere <- sum(ari_MCEG_hemisphere > ari_ZG3_hemisphere))
(pvalue_MCEG_ZG3_hemisphere <- 1 - pbinom(win_MCEG_ZG3_hemisphere, ngraphs, 0.5))
(win_MCEG_Louvain_hemisphere <- sum(ari_MCEG_hemisphere > ari_Louvain_hemisphere))
(pvalue_MCEG_Louvain_hemisphere <- 1 - pbinom(win_MCEG_Louvain_hemisphere, ngraphs, 0.5))
(win_MCEG_Rw_hemisphere <- sum(ari_MCEG_hemisphere > ari_Rw_hemisphere))
(pvalue_MCEG_Rw_hemisphere <- 1 - pbinom(win_MCEG_Rw_hemisphere, ngraphs, 0.5))

(win_MCG_ZG1_tissue <- sum(ari_MCG_tissue > ari_ZG1_tissue))
(pvalue_MCG_ZG1_tissue <- 1 - pbinom(win_MCG_ZG1_tissue, ngraphs, 0.5))
(win_MCG_ZG2_tissue <- sum(ari_MCG_tissue > ari_ZG2_tissue))
(pvalue_MCG_ZG2_tissue <- 1 - pbinom(win_MCG_ZG2_tissue, ngraphs, 0.5))
(win_MCG_ZG3_tissue <- sum(ari_MCG_tissue > ari_ZG3_tissue))
(pvalue_MCG_ZG3_tissue <- 1 - pbinom(win_MCG_ZG3_tissue, ngraphs, 0.5))
(win_MCG_Louvain_tissue <- sum(ari_MCG_tissue > ari_Louvain_tissue))
(pvalue_MCG_Louvain_tissue <- 1 - pbinom(win_MCG_Louvain_tissue, ngraphs, 0.5))
(win_MCG_Rw_tissue <- sum(ari_MCG_tissue > ari_Rw_tissue))
(pvalue_MCG_Rw_tissue <- 1 - pbinom(win_MCG_Rw_tissue, ngraphs, 0.5))

(win_MCEG_ZG1_tissue <- sum(ari_MCEG_tissue > ari_ZG1_tissue))
(pvalue_MCEG_ZG1_tissue <- 1 - pbinom(win_MCEG_ZG1_tissue, ngraphs, 0.5))
(win_MCEG_ZG2_tissue <- sum(ari_MCEG_tissue > ari_ZG2_tissue))
(pvalue_MCEG_ZG2_tissue <- 1 - pbinom(win_MCEG_ZG2_tissue, ngraphs, 0.5))
(win_MCEG_ZG3_tissue <- sum(ari_MCEG_tissue > ari_ZG3_tissue))
(pvalue_MCEG_ZG3_tissue <- 1 - pbinom(win_MCEG_ZG3_tissue, ngraphs, 0.5))
(win_MCEG_Louvain_tissue <- sum(ari_MCEG_tissue > ari_Louvain_tissue))
(pvalue_MCEG_Louvain_tissue <- 1 - pbinom(win_MCEG_Louvain_tissue, ngraphs, 0.5))
(win_MCEG_Rw_tissue <- sum(ari_MCEG_tissue > ari_Rw_tissue))
(pvalue_MCEG_Rw_tissue <- 1 - pbinom(win_MCEG_Rw_tissue, ngraphs, 0.5))

(win_MCG_ZG1_Y <- sum(ari_MCG_Y > ari_ZG1_Y))
(pvalue_MCG_ZG1_Y <- 1 - pbinom(win_MCG_ZG1_Y, ngraphs, 0.5))
(win_MCG_ZG2_Y <- sum(ari_MCG_Y > ari_ZG2_Y))
(pvalue_MCG_ZG2_Y <- 1 - pbinom(win_MCG_ZG2_Y, ngraphs, 0.5))
(win_MCG_ZG3_Y <- sum(ari_MCG_Y > ari_ZG3_Y))
(pvalue_MCG_ZG3_Y <- 1 - pbinom(win_MCG_ZG3_Y, ngraphs, 0.5))
(win_MCG_Louvain_Y <- sum(ari_MCG_Y > ari_Louvain_Y))
(pvalue_MCG_Louvain_Y <- 1 - pbinom(win_MCG_Louvain_Y, ngraphs, 0.5))
(win_MCG_Rw_Y <- sum(ari_MCG_Y > ari_Rw_Y))
(pvalue_MCG_Rw_Y <- 1 - pbinom(win_MCG_Rw_Y, ngraphs, 0.5))

(win_MCEG_ZG1_Y <- sum(ari_MCEG_Y > ari_ZG1_Y))
(pvalue_MCEG_ZG1_Y <- 1 - pbinom(win_MCEG_ZG1_Y, ngraphs, 0.5))
(win_MCEG_ZG2_Y <- sum(ari_MCEG_Y > ari_ZG2_Y))
(pvalue_MCEG_ZG2_Y <- 1 - pbinom(win_MCEG_ZG2_Y, ngraphs, 0.5))
(win_MCEG_ZG3_Y <- sum(ari_MCEG_Y > ari_ZG3_Y))
(pvalue_MCEG_ZG3_Y <- 1 - pbinom(win_MCEG_ZG3_Y, ngraphs, 0.5))
(win_MCEG_Louvain_Y <- sum(ari_MCEG_Y > ari_Louvain_Y))
(pvalue_MCEG_Louvain_Y <- 1 - pbinom(win_MCEG_Louvain_Y, ngraphs, 0.5))
(win_MCEG_Rw_Y <- sum(ari_MCEG_Y > ari_Rw_Y))
(pvalue_MCEG_Rw_Y <- 1 - pbinom(win_MCEG_Rw_Y, ngraphs, 0.5))

# print(load("RData/df-all-drosophila.RData"))
# df.all <- df.all %>%
#     mutate(SMS_tissue = ari_MCG_tissue, SMS_hemisphere = ari_MCG_hemisphere, SMS_Y = ari_MCG_Y, SMS_region = ari_MCG_region,
#            SMS_R_tissue = ari_MCEG_tissue, SMS_R_hemisphere = ari_MCEG_hemisphere, SMS_R_Y = ari_MCEG_Y, SMS_R_region = ari_MCEG_region,
#            ZG1_tissue = ari_ZG1_tissue, ZG1_hemisphere = ari_ZG1_hemisphere, ZG1_Y = ari_ZG1_Y, ZG1_region = ari_ZG1_region,
#            ZG2_tissue = ari_ZG2_tissue, ZG2_hemisphere = ari_ZG2_hemisphere, ZG2_Y = ari_ZG2_Y, ZG2_region = ari_ZG2_region,
#            ZG3_tissue = ari_ZG3_tissue, ZG3_hemisphere = ari_ZG3_hemisphere, ZG3_Y = ari_ZG3_Y, ZG3_region = ari_ZG3_region)
#save(df.all, file="RData/df-all-drosophila.RData")
print(load("RData/df-all-drosophila.RData"))

(win_MCG_IRM_hemisphere <- sum(ari_MCG_hemisphere > df.all$IRM_hemispehre))
(pvalue_MCG_IRM_hemisphere <- 1 - pbinom(win_MCG_IRM_hemisphere, ngraphs, 0.5))
(win_MCEG_IRM_hemisphere <- sum(ari_MCEG_hemisphere > df.all$IRM_hemispehre))
(pvalue_MCEG_IRM_hemisphere <- 1 - pbinom(win_MCEG_IRM_hemisphere, ngraphs, 0.5))

(win_MCG_IRM_tissue <- sum(ari_MCG_tissue > df.all$IRM_tissue))
(pvalue_MCG_IRM_tissue <- 1 - pbinom(win_MCG_IRM_tissue, ngraphs, 0.5))
(win_MCEG_IRM_tissue <- sum(ari_MCEG_tissue > df.all$IRM_hemispehre))
(pvalue_MCEG_IRM_tissue <- 1 - pbinom(win_MCEG_IRM_tissue, ngraphs, 0.5))

(win_MCG_IRM_Y <- sum(ari_MCG_Y > df.all$IRM_Y))
(pvalue_MCG_IRM_Y <- 1 - pbinom(win_MCG_IRM_Y, ngraphs, 0.5))
(win_MCEG_IRM_Y <- sum(ari_MCEG_Y > df.all$IRM_hemispehre))
(pvalue_MCEG_IRM_Y <- 1 - pbinom(win_MCEG_IRM_Y, ngraphs, 0.5))

