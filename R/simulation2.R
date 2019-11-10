source("loadData.R")

requiredPkg <- c("googledrive", "tidyverse")
newPkg <- requiredPkg[!(requiredPkg %in% installed.packages()[, "Package"])]
if (length(newPkg)) install.packages(newPkg, dependencies = TRUE)
sapply(requiredPkg, require, character.only = TRUE)

findex <- c("09", "095", "1", "105", "11", "115")
meanAriB2 <- rep(0,6)
meanAriB2s <- rep(0,6)
meanAriSVT1 <- rep(0,6)
meanAriSVT2 <- rep(0,6)
meanAriSVT3 <- rep(0,6)
meanDB2s <- rep(0,6)
meanDSVT1 <- rep(0,6)
meanDSVT2 <- rep(0,6)
meanDSVT3 <- rep(0,6)

drive_auth_config(active = FALSE)
folder_url <- "https://drive.google.com/open?id=1PXyy6-z_O-ICu79kjZB-W9hkr-vgN6kh"
folder <- drive_get(as_id(folder_url))
files <- drive_ls(folder)

##############################################
# plot figure 8
##############################################
files1 <- files %>% filter(grepl(glob2rx("eig_BF1_n500_p0095*"), name))
result <- loadData(files1)
# hist(result$ariBlock2sList -result$ariSvt1List, breaks=seq(-0.05, 0.2, 0.02), prob=T,col=1,density=10,angle=135,
#      xlim = c(-0.05,0.2), ylim = c(0, 40), main = "", xlab = "Diffence of ARI", lty = 1)
# lines(stats::density(result$ariBlock2sList -result$ariSvt1List,adj=2),col=1,lwd=2)
# hist(result$ariBlock2sList -result$ariSvt2List, breaks=seq(-0.05, 0.2, 0.02), prob=T,col=2,density=10,angle=90,
#      add = TRUE, lty = 1)
# lines(stats::density(result$ariBlock2sList -result$ariSvt2List,adj=10),col=2,lwd=2)
# hist(result$ariBlock2sList -result$ariSvt3List, breaks=seq(-0.05, 0.2, 0.02), prob=T,col=4,density=10,angle=45,
#      add = TRUE, lty = 1)
# lines(stats::density(result$ariBlock2sList -result$ariSvt3List,adj=2),col=4,lwd=2)
# legend("topleft", col = c(1,2,4), legend = c("MCG-ZG1", "MCG-ZG2", "MCG-ZG3"),
#        lwd = c(1.5, 1.5, 1.5))

df.a <- data.frame(p="p = 0.095",
                   ZG1=result$ariBlock2sList - result$ariSvt1List,
                   ZG2=result$ariBlock2sList - result$ariSvt2List,
                   ZG3=result$ariBlock2sList - result$ariSvt3List)
df.a <- df.a %>% gather(`ZG1`,`ZG2`,`ZG3`, key="method", value=diff)
# df.a %>% ggplot(aes(x=method,y=diff,color=method,fill=method)) +
#   geom_boxplot(alpha=0.5,notch = TRUE) +
#   theme(legend.position="none") +
#   ylab("Difference of ARI")

##############################################
files2 <- files %>% filter(grepl(glob2rx("eig_BF1_n500_p0115*"), name))
result <- loadData(files2)
# hist(result$ariBlock2sList -result$ariSvt1List, breaks=seq(-0.09, 0.19, 0.02), prob=T,col=1,density=10,angle=135,
#      xlim = c(-0.1,0.2), ylim = c(0, 20), main = "", xlab = "Diffence of ARI", lty = 1)
# lines(stats::density(result$ariBlock2sList -result$ariSvt1List,adj=2),col=1,lwd=2)
# hist(result$ariBlock2sList -result$ariSvt2List, breaks=seq(-0.09, 0.19, 0.02), prob=T,col=2,density=10,angle=90,
#      add = TRUE, lty = 1)
# lines(stats::density(result$ariBlock2sList -result$ariSvt2List,adj=2),col=2,lwd=2)
# hist(result$ariBlock2sList -result$ariSvt3List, breaks=seq(-0.09, 0.19, 0.02), prob=T,col=4,density=10,angle=45,
#      add = TRUE, lty = 1)
# lines(stats::density(result$ariBlock2sList -result$ariSvt3List,adj=2),col=4,lwd=2)
# legend("topleft", col = c(1,2,4), legend = c("MCG-ZG1", "MCG-ZG2", "MCG-ZG3"),
#        lwd = c(1.5, 1.5, 1.5))

df.b <- data.frame(p="p = 0.115",
                   ZG1=result$ariBlock2sList - result$ariSvt1List,
                   ZG2=result$ariBlock2sList - result$ariSvt2List,
                   ZG3=result$ariBlock2sList - result$ariSvt3List)
df.b <- df.b %>% gather(`ZG1`,`ZG2`,`ZG3`, key="method", value=diff)
# df.b %>% ggplot(aes(x=method,y=diff,color=method,fill=method)) +
#   geom_boxplot(alpha=0.5,notch = TRUE) +
#   theme(legend.position="none") +
#   ylab("Difference of ARI")

df.ab <- rbind(df.a, df.b)
#df.ab$p <- factor(df.ab$p)
df.ab %>% ggplot(aes(x=method,y=diff,color=method,fill=method)) +
  facet_grid(~p) +
  geom_boxplot(alpha=0.5,notch = TRUE) +
#  geom_violin(alpha=0.5) +
  theme(legend.position="none") +
  ylab("ARI(MCEG) - ARI(ZG)")
#ggsave("Fig8-alt.pdf", width=9, height=5)

##############################################
for (i in 1:6){
  filename = paste("eig_BF1_n500_p0", findex[i], "_mc*", sep="")
  filei <- files %>% filter(grepl(glob2rx(filename), name))

  result <- loadData(filei)

  meanAriB2[i] <- mean(result$ariBlock2List)
  meanAriB2s[i] <- mean(result$ariBlock2sList)
  meanAriSVT1[i] <- mean(result$ariSvt1List)
  meanAriSVT2[i] <- mean(result$ariSvt2List)
  meanAriSVT3[i] <- mean(result$ariSvt3List)

  meanDB2s[i] <- mean(result$dhatBlock2List)
  meanDSVT1[i] <- mean(result$dhatSvt1List)
  meanDSVT2[i] <- mean(result$dhatSvt2List)
  meanDSVT3[i] <- mean(result$dhatSvt3List)
}

##############################################
# plot figure 10
xrange <- c(0.09,0.115)
yrange <- c(0,1)
plot(xrange, yrange, type = "n", xlab = "p", ylab = "ARI", col.axis = "black")
x <- c(0.09, 0.095,0.1,0.105,0.11,0.115)
y <- rep(0,100)
# the results of Louvain and Walktrap are obtained separately
meanAriLouvain <- c(0.094, 0.052, 0.033, 0.018, 0.013, 0.008)
meanAriWalktrap <- c(0.787, 0.633, 0.381, 0.203, 0.115, 0.070)

lines(x, meanAriSVT1, type = "l", lwd = 1.5, lty = 1, col = 1)
lines(x, meanAriSVT2, type = "l", lwd = 1.5, lty = 1, col = 4)
lines(x, meanAriSVT3, type = "l", lwd = 1.5, lty = 1, col = 2)
lines(x, meanAriB2, type = "l", lwd = 1.5, lty = 1, col = 8)
lines(x, meanAriB2s, type = "l", lwd = 1.5, lty = 1, col = 6)
lines(x, meanAriLouvain, type = "l", lwd = 1.5, lty = 1, col = 3)
lines(x, meanAriWalktrap, type = "l", lwd = 1.5, lty = 1, col = 5)

legend("bottomright", col = c(1,4,2,8,6,3,5),
       legend = c("ZG1", "ZG2", "ZG3", "MCG", "MCEG", "Louvain", "Walktrap"),
       lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5), cex = 0.5)
#dev.print(pdf, "simulation_all_methods.pdf", width=5, height=4)
##############################################
# plot figure 11(a)
op <- par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

xrange <- c(0.09,0.115)
yrange <- c(0.8,1)
plot(xrange, yrange, type = "n", xlab = "p", ylab = "mean(ARI)", col.axis = "black")
x <- c(0.09, 0.095,0.1,0.105,0.11,0.115)
y <- rep(0,100)

lines(x, meanAriSVT1, type = "l", lwd = 3, lty = 1, col = 1)
lines(x, meanAriSVT2, type = "l", lwd = 3, lty = 1, col = 4)
lines(x, meanAriSVT3, type = "l", lwd = 3, lty = 1, col = 2)
lines(x, meanAriB2s, type = "l", lwd = 4.5, lty = 1, col = 6)

legend("topright", col = c(1,4,2,6),
       legend = c("ZG1", "ZG2", "ZG3", "MCEG"),
       lwd = c(1.5, 1.5, 1.5, 1.5), cex = 0.7)
#dev.print(pdf, "simulation3_arimean.pdf")

##############################################
# plot figure 11(b)
xrange <- c(0.09,0.115)
yrange <- c(-1.5,4)
plot(xrange, yrange, type = "n", xlab = "p", ylab = expression(paste("mean(", hat(d), " - ", d, ")")), col.axis = "black")
x <- c(0.09, 0.095,0.1,0.105,0.11,0.115)
y <- rep(0,100)

lines(x, meanDSVT1-2, type = "l", lwd = 3, lty = 1, col = 1)
lines(x, meanDSVT2-2, type = "l", lwd = 3, lty = 1, col = 4)
lines(x, meanDSVT3-2, type = "l", lwd = 3, lty = 1, col = 2)
lines(x, meanDB2s-2, type = "l", lwd = 4.5, lty = 1, col = 6)
lines(x, rep(0,6), type = "l", lwd = 1.5, lty = 2, col = 3)

legend("topright", col = c(1,4,2,6),
       legend = c("ZG1", "ZG2", "ZG3", "MCEG"),
       lwd = c(1.5, 1.5, 1.5, 1.5), cex = 0.7)
#dev.print(pdf, "simulation3_Dmean.pdf")
par(op)
