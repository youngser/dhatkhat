source("loadData.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

requiredPkg <- c("XML", "tidyverse")
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

# drive_auth_config(active = FALSE)
# folder_url <- "https://drive.google.com/open?id=1PXyy6-z_O-ICu79kjZB-W9hkr-vgN6kh"
# folder <- drive_get(as_id(folder_url))
# files <- drive_ls(folder)

folder_url <- "http://www.cis.jhu.edu/~parky/dhatKhat/Results/result_simulation/"
docs <- XML::htmlParse(folder_url)
links <- XML::xpathSApply(docs, "//a/@href")
dfiles <- links[grepl(glob2rx("*.RData"), links)]
files <- paste0(folder_url, dfiles)
#print(load(url(files[x])))


##############################################
# plot figure 8
##############################################

#files1 <- files %>% filter(grepl(glob2rx("eig_BF1_n500_p0095*"), name))
files1 <- files[grepl(glob2rx("eig_BF1_n500_p0095*"),dfiles)]
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

df8.a <- data.frame(p="p = 0.095",
                    SMS = result$ariBlock2sList,
                    ZG1 = result$ariSvt1List,
                    ZG2 = result$ariSvt2List,
                    ZG3 = result$ariSvt3List) #%>%
#        gather(key = "group", value = "ari", -p)


##############################################
#files2 <- files %>% filter(grepl(glob2rx("eig_BF1_n500_p0115*"), name))
files2 <- files[grepl(glob2rx("eig_BF1_n500_p0115*"),dfiles)]
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
df8.b <- data.frame(p="p = 0.115",
                    SMS = result$ariBlock2sList,
                    ZG1 = result$ariSvt1List,
                    ZG2 = result$ariSvt2List,
                    ZG3 = result$ariSvt3List) #%>%
#        gather(key = "group", value = "ari", -p)

df.ab <- rbind(df.a, df.b)
#df.ab$p <- factor(df.ab$p)
df.ab %>% ggplot(aes(x=method,y=diff,color=method,fill=method)) +
  facet_grid(~p) +
  geom_boxplot(alpha=0.5,notch = TRUE) +
#  geom_violin(alpha=0.5) +
  theme(legend.position="none") +
  ylab("ARI(SMS) - ARI(ZG)") +
  theme(axis.title=element_text(size=15)) +
  theme(axis.text.x=element_text(size=15)) +
  #  theme(axis.ticks = element_line(size = 5)) +
  theme(axis.text.y=element_text(size=15)) +
  theme(strip.text=element_text(size=rel(1.2))) +
  theme(legend.text = element_text(colour="black", size = 15, face = "plain"))

ggsave("Fig8-alt.pdf", width=9, height=5)

## paired t.test
t.test(as.numeric(df8.a$SMS), as.numeric(df8.a$ZG1), paired=TRUE, alternative = "gre")
t.test(as.numeric(df8.a$SMS), as.numeric(df8.a$ZG2), paired=TRUE, alternative = "gre")
t.test(as.numeric(df8.a$SMS), as.numeric(df8.a$ZG3), paired=TRUE, alternative = "gre")

t.test(as.numeric(df8.b$SMS), as.numeric(df8.b$ZG1), paired=TRUE, alternative = "gre")
t.test(as.numeric(df8.b$SMS), as.numeric(df8.b$ZG2), paired=TRUE, alternative = "gre")
t.test(as.numeric(df8.b$SMS), as.numeric(df8.b$ZG3), paired=TRUE, alternative = "gre")


##############################################
# plot figure 10
##############################################

for (i in 1:6){
  filename = paste("eig_BF1_n500_p0", findex[i], "_mc*", sep="")
#  filei <- files %>% filter(grepl(glob2rx(filename), name))
  filei <- files[grepl(glob2rx(filename),dfiles)]
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

x <- c(0.09, 0.095,0.1,0.105,0.11,0.115)
y <- rep(0,100)
# the results of Louvain and Walktrap are obtained separately, calculated from `LouvainAndRw.R`
meanAriLouvain <- c(0.094, 0.052, 0.033, 0.018, 0.013, 0.008)
meanAriWalktrap <- c(0.787, 0.633, 0.381, 0.203, 0.115, 0.070)
meanAriIRM <- c(0.132, 0.128, 0.125, 0.121, 0.116, 0.111)

mycol <- gg_color_hue(8)

xrange <- c(0.09,0.115)
yrange <- c(0,1)

plot(xrange, yrange, type = "n", xlab = "p", ylab = "ARI", col.axis = "black")
lines(x, meanAriSVT1, type = "l", lwd = 1.5, lty = 1, col = mycol[1]) # ZG1
lines(x, meanAriSVT2, type = "l", lwd = 1.5, lty = 1, col = mycol[2]) # ZG2
lines(x, meanAriSVT3, type = "l", lwd = 1.5, lty = 1, col = mycol[3]) # ZG3
lines(x, meanAriB2, type = "l", lwd = 1.5, lty = 1, col = mycol[4])   # SMS
lines(x, meanAriB2s, type = "l", lwd = 1.5, lty = 1, col = mycol[5])  # SMS_Reduced
lines(x, meanAriLouvain, type = "l", lwd = 1.5, lty = 1, col = mycol[6])
lines(x, meanAriWalktrap, type = "l", lwd = 1.5, lty = 1, col = mycol[7])
lines(x, meanAriIRM, type = "l", lwd = 1.5, lty = 1, col = mycol[8])
legend(0.108, 0.6, col = mycol,
       legend = c("ZG1", "ZG2", "ZG3", "SMS", "SMS_Reduced", "Louvain", "Walktrap", "IRM"),
       lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5), cex = 0.5, y.intersp=1.5)
box()
dev.print(pdf, "simulation_all_methods.pdf", width=5, height=4)


##############################################
# plot figure 11(a)
##############################################

op <- par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

xrange <- c(0.09,0.115)
yrange <- c(0.8,1)
plot(xrange, yrange, type = "n", xlab = "p", ylab = "mean(ARI)", col.axis = "black")
x <- c(0.09, 0.095,0.1,0.105,0.11,0.115)
y <- rep(0,100)

lines(x, meanAriSVT1, type = "l", lwd = 3, lty = 1, col = mycol[1])
lines(x, meanAriSVT2, type = "l", lwd = 3, lty = 1, col = mycol[2])
lines(x, meanAriSVT3, type = "l", lwd = 3, lty = 1, col = mycol[3])
lines(x, meanAriB2s, type = "l", lwd = 4.5, lty = 1, col = mycol[5])

legend("topright", col = mycol[c(1:3,5)],#c(1,4,2,6),
       legend = c("ZG1", "ZG2", "ZG3", "SMS_Reduced"),
       lwd = c(1.5, 1.5, 1.5, 1.5), cex = 0.7)
dev.print(pdf, "simulation3_arimean.pdf")

##############################################
# plot figure 11(b)
##############################################

xrange <- c(0.09,0.115)
yrange <- c(-1.5,4)
plot(xrange, yrange, type = "n", xlab = "p", ylab = expression(paste("mean(", hat(d), " - ", d, ")")), col.axis = "black")
x <- c(0.09, 0.095,0.1,0.105,0.11,0.115)
y <- rep(0,100)

lines(x, meanDSVT1-2, type = "l", lwd = 3, lty = 1, col = mycol[1])
lines(x, meanDSVT2-2, type = "l", lwd = 3, lty = 1, col = mycol[2])
lines(x, meanDSVT3-2, type = "l", lwd = 3, lty = 1, col = mycol[3])
lines(x, meanDB2s-2, type = "l", lwd = 4.5, lty = 1, col = mycol[5])
lines(x, rep(0,6), type = "l", lwd = 1.5, lty = 2, col = 1)

# legend("topright", col = mycol[c(1:3,5)],#c(1,4,2,6),
#        legend = c("ZG1", "ZG2", "ZG3", "SMS_Reduced"),
#        lwd = c(1.5, 1.5, 1.5, 1.5), cex = 0.7)
dev.print(pdf, "simulation3_Dmean.pdf")
par(op)
