source("ase.R")
source("emase.R")
source("clt.R")

requiredPkg <- c("googledrive", "tidyverse", "parallel", "ellipse")
newPkg <- requiredPkg[!(requiredPkg %in% installed.packages()[, "Package"])]
if (length(newPkg)) install.packages(newPkg, dependencies = TRUE)
sapply(requiredPkg, require, character.only = TRUE)

drive_auth_config(active = FALSE)
folder_url <- "https://drive.google.com/open?id=1SD2RLcgww_rvr53l0KDQQBe38RncvULN"
folder <- drive_get(as_id(folder_url))
files <- drive_ls(folder)

# call "figure1(K=2)" for Fig 1 (a),
# call "figure1(K=3)" for Fig 1 (b)
figure1 <- function(n = 1000, K = 2){
    set.seed(1240)

    if (K == 2) {
        ## inf2.pdf
        rho <- c(0.5,0.5)
        B <- matrix(c(0.2,0.1,0.1,0.25), nrow = 2, ncol = 2)
    } else {
        ## inf3.pdf
        rho <- c(0.4,0.4,0.2)
        B <- rbind(c(0.2,0.1,0.08),c(0.1,0.25,0.05), c(0.08,0.05,0.4))
    }
    #    K <- length(rho)
    nB <- n*rho
    tau <- rep(1:K, times=nB)

    #    A <- rg.SBM(n, B, rho)
    g <- sample_sbm(n, B, nB, directed=FALSE)
    #    tau <- A$tau
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))
    X <- xx[tau,]

    #    Xhat <- stfp(A$adjacency, dim=K)
    Xhat <- ase(g, max.d=K, scale.by.values=TRUE)
    W <- procrustes(Xhat, X)
    residue <- sqrt(n)*(Xhat%*%W - X)
    print(norm(Xhat%*%W - X,"F")^2)

    sigs <- block.var(xx, rho)
    print(sigs)

    df <- data.frame(x = (Xhat%*%W)[,1], y = (Xhat%*%W)[,2],
                     block = as.character(tau),
                     tau = tau)

    df_ell <- data.frame()
    for(g in unique(df$tau)){
        sig <- sigs[[g]]
        df_ell <- rbind(df_ell, cbind(as.data.frame(
            with(df[df$tau==g,],
                 ellipse(sig[1,2]/sqrt(sig[1,1]*sig[2,2]),
                         scale=c(sqrt(sig[1,1])/sqrt(n),sqrt(sig[2,2])/sqrt(n)),
                         centre=c(xx[g,1],xx[g,2]))),group=g)))

    }

    df_ell <- cbind(df_ell, block=as.character(rep(1:K, each=nrow(df_ell)/K)))
    df.xx <- as.data.frame(xx[,1:2]); names(df.xx) <- c("x","y")
    df.xx <- cbind(df.xx, block=as.character(c(1:K)))

    ggplot(data=df, aes(x=x, y=y, color=block)) +
        geom_point(size=2, alpha=0.5) +
        stat_ellipse(type="norm", size=1, linetype = "dashed", aes(color=block)) +
        geom_point(data=df.xx, aes(x=x,y=y), color="black", size=3) +
        geom_path(data=df_ell, aes(x=x, y=y, group=block), size=1, linetype = 2, colour=1) +
        #    geom_path(data = df_ell[-c(1:(nrow(df_ell)/K)),], aes(x = x, y = y), size = 1, linetype = 2, colour = 1) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
        #  theme_light() +
        theme(legend.position = "none")

    #    ggsave(paste0("inf",K,".pdf"), width=5, height=5)
}

figure23data <- function()
{
    N <- 100
    seed <- 0
    ns <- c(200,500,1000,2000,4000,8000,16000)
    B2 <- rbind(c(0.2,0.1),c(0.1,0.25))
    B3 <- rbind(c(0.2,0.1,0.08),c(0.1,0.25,0.05),
                c(0.08,0.05,0.4))
    pi2 <- c(0.5,0.5)
    pi3 <- c(0.4,0.4,0.2)
    mus2 <- vector('list',length(ns))
    mus3 <- vector('list',length(ns))
    names(mus2) <- ns
    names(mus3) <- ns
    vars2 <- vector('list',length(ns))
    vars3 <- vector('list',length(ns))
    names(vars2) <- ns
    names(vars3) <- ns
    varsS2 <- vector('list',length(ns))
    varsS3 <- vector('list',length(ns))
    names(varsS2) <- ns
    names(varsS3) <- ns
    V <- vector('list',length(ns))
    names(V) <- ns
    VS <- vector('list',length(ns))
    names(VS) <- ns
    for(i in 1:length(ns)){
        n <- ns[i]
        cat(n,"\n")
        ng2 <- n*pi2
        ng3 <- n*pi3
        cls2 <- rep(1:2,times=ng2)
        cls3 <- rep(1:3,times=ng3)
        mus2[[i]] <- vector('list',2)
        mus2[[i]][[1]] <- matrix(NA,nrow=N,ncol=80)
        mus2[[i]][[2]] <- matrix(NA,nrow=N,ncol=80)
        mus3[[i]] <- vector('list',3)
        mus3[[i]][[1]] <- matrix(NA,nrow=N,ncol=80)
        mus3[[i]][[2]] <- matrix(NA,nrow=N,ncol=80)
        mus3[[i]][[3]] <- matrix(NA,nrow=N,ncol=80)
        vars2[[i]] <- vector('list',2)
        vars2[[i]][[1]] <- matrix(NA,nrow=N,ncol=80)
        vars2[[i]][[2]] <- matrix(NA,nrow=N,ncol=80)
        V[[i]] <- matrix(NA,nrow=N,ncol=80)
        VS[[i]] <- matrix(NA,nrow=N,ncol=80)
        vars3[[i]] <- vector('list',3)
        vars3[[i]][[1]] <- matrix(NA,nrow=N,ncol=80)
        vars3[[i]][[2]] <- matrix(NA,nrow=N,ncol=80)
        vars3[[i]][[3]] <- matrix(NA,nrow=N,ncol=80)
        varsS2[[i]] <- vector('list',2)
        varsS2[[i]][[1]] <- matrix(NA,nrow=N,ncol=80)
        varsS2[[i]][[2]] <- matrix(NA,nrow=N,ncol=80)
        varsS3[[i]] <- vector('list',3)
        varsS3[[i]][[1]] <- matrix(NA,nrow=N,ncol=80)
        varsS3[[i]][[2]] <- matrix(NA,nrow=N,ncol=80)
        varsS3[[i]][[3]] <- matrix(NA,nrow=N,ncol=80)
        for(j in 1:N){
            set.seed(seed+i*N*2+j)
            g2 <- sample_sbm(n,B2,ng2,directed=FALSE)
            emb2 <- ase(g2,verbose=FALSE,max.d=80,scale.by.values=FALSE)
            embS2 <- ase(g2,verbose=FALSE,max.d=80,scale.by.values=TRUE)
            mus2[[i]][[1]][j,] <- apply(emb2[which(cls2==1),],2,mean)
            mus2[[i]][[2]][j,] <- apply(emb2[which(cls2==2),],2,mean)
            vars2[[i]][[1]][j,] <- apply(emb2[which(cls2==1),],2,var)
            vars2[[i]][[2]][j,] <- apply(emb2[which(cls2==2),],2,var)
            varsS2[[i]][[1]][j,] <- apply(embS2[which(cls2==1),],2,var)
            varsS2[[i]][[2]][j,] <- apply(embS2[which(cls2==2),],2,var)
            V[[i]][j,] <- apply(emb2,2,var)
            VS[[i]][j,] <- apply(embS2,2,var)
            g3 <- sample_sbm(n,B3,ng3,directed=FALSE)
            emb3 <- ase(g3,verbose=FALSE,max.d=80,scale.by.values=FALSE)
            embS3 <- ase(g3,verbose=FALSE,max.d=80,scale.by.values=TRUE)
            mus3[[i]][[1]][j,] <- apply(emb3[which(cls3==1),],2,mean)
            mus3[[i]][[2]][j,] <- apply(emb3[which(cls3==2),],2,mean)
            mus3[[i]][[3]][j,] <- apply(emb3[which(cls3==3),],2,mean)
            vars3[[i]][[1]][j,] <- apply(emb3[which(cls3==1),],2,var)
            vars3[[i]][[2]][j,] <- apply(emb3[which(cls3==2),],2,var)
            vars3[[i]][[3]][j,] <- apply(emb3[which(cls3==3),],2,var)
            varsS3[[i]][[1]][j,] <- apply(embS3[which(cls3==1),],2,var)
            varsS3[[i]][[2]][j,] <- apply(embS3[which(cls3==2),],2,var)
            varsS3[[i]][[3]][j,] <- apply(embS3[which(cls3==3),],2,var)
        }
        save(i,n,mus2,mus3,vars2,vars3,V,varsS2,varsS3,VS,N,ns,seed,
             B2,B3,pi2,pi3,
             file=paste0('fig2_3',n,'.RData'))
    }

    save(mus2,mus3,vars2,vars3,V,varsS2,varsS3,VS,N,ns,seed,
         B2,B3,pi2,pi3,
         file='fig2_3.RData')
}

# Fig 2
figure2 <- function(max.d=80)
{
    files1 <- files %>% filter(grepl(glob2rx("fig2_3*"), name))
    file <- drive_download(files1, overwrite = TRUE, verbose = FALSE)
    load(file$local_path)
#    load('fig2_3.RData')

    max.d <- min(max.d,ncol(mus2[[1]][[1]]))

    m <- par('mar')
    par(mar=c(m[1:2],1,1))
    par(mfrow=c(2,2))
    boxplot(mus2[[1]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[1]), cex.lab=1.5)
    axis(2)
    axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
    title("n=200")
    box()
    boxplot(mus2[[4]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[1]), cex.lab=1.5)
    axis(2)
    axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
    title("n=2000")
    box()

    boxplot(mus2[[1]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[2]), cex.lab=1.5)
    axis(2)
    axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
    title("n=200")
    box()
    boxplot(mus2[[4]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[2]), cex.lab=1.5)
    axis(2)
    axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
    title("n=2000")
    box()
#    dev.print(device=pdf,file='fig1means.pdf', width=8, height=8) # Fig 2 (a)
    par(mar=m)

    m <- par('mar')
    par(mfrow=c(2,1))
    par(mar=c(m[1:2],0.5,1))
    M1 <- lapply(mus2,function(x) (c(x[[1]][,-(1:2)])))
    M2 <- lapply(mus2,function(x) (c(x[[2]][,-(1:2)])))
    boxplot(M1,names=ns,notch=TRUE,xlab="n",ylab=expression(mu[1]), cex.lab=1.5)
    boxplot(M2,names=ns,notch=TRUE,xlab="n",ylab=expression(mu[2]), cex.lab=1.5)
    par(mar=m)
#    dev.print(device=pdf,file='musall.pdf', width=8, height=8) # Fig 2 (b)
    par(mfrow=c(1,1))

}

figure3 <- function(max.d=80,useS=FALSE)
{
    files1 <- files %>% filter(grepl(glob2rx("fig2_3*"), name))
    file <- drive_download(files1, overwrite = TRUE, verbose = FALSE)
    load(file$local_path)

    max.d <- min(max.d,ncol(mus2[[1]][[1]]))
    # if(useS){
    #     vars2 <- varsS2
    #     vars3 <- varsS3
    #     V <- VS
    # }

    par(mfrow=c(2,1))
    m <- par('mar')
    par(mar=c(m[1],m[2]+1,0.75,m[4]))
    ylim <- range(unlist(vars2[[3]]))
    boxplot(vars2[[3]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,
            ylim=ylim,
            xlab='Dimension',ylab=expression(sigma^2), cex.lab=1.5)
    axis(2)
    axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
    #title("n=200")
    box()
    boxplot(vars2[[3]][[2]][,-(1:2)],add=TRUE,names=3:max.d,axes=FALSE, cex.lab=1.5)

    ylim <- range(unlist(vars2[[7]]))
    boxplot(vars2[[7]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,
            ylim=ylim,
            xlab='Dimension',ylab=expression(sigma^2), cex.lab=1.5)
    axis(2)
    axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
    #title("n=16000")
    box()
    boxplot(vars2[[7]][[2]][,-(1:2)],add=TRUE,names=3:max.d,axes=FALSE, cex.lab=1.5)
#    dev.print(device=pdf,file='fig2sigmas16000.pdf') # Fig 3 (a)
    par(mar=m)

    par(mfrow=c(1,1))
    yl <- 1.35
    cols <- 1:length(vars2)
    for(i in 1:length(vars2)){
        yl <- range(yl,ns[i]*apply(vars2[[i]][[1]][,-(1:2)],2,median))
        yl <- range(yl,ns[i]*apply(vars2[[i]][[2]][,-(1:2)],2,median))
    }
    plot(3:80,ylim=yl,xlab="Dimension",
         ylab=expression(n*sigma^2),type='n',
         xlim=c(3,82), cex.lab=1.5)
    for(i in 1:length(vars2)){
        data <- ns[i]*apply(vars2[[i]][[1]][,-(1:2)],2,median)
        lines(3:80,data,lty=1,col=cols[i],lwd=2)
        data <- ns[i]*apply(vars2[[i]][[2]][,-(1:2)],2,median)
        lines(3:80,data,lty=2,col=cols[i],lwd=2)
    }
#    legend(0,yl[2]+.025,legend=paste0("n=",ns),ncol=4,col=cols,lty=1,lwd=1)
#    legend(70,1.05,legend=c("B1","B2"),lty=1:2,col=1)
#    dev.print(device=pdf,file='allsigmas.pdf') # Fig 3 (b)
}

figure5Data <- function()
{
    N <- 100
    seed <- 0
    ns <- c(200,500,1000,2000,4000)
    B2 <- rbind(c(0.2,0.1),c(0.1,0.25))
    B3 <- rbind(c(0.2,0.1,0.08),c(0.1,0.25,0.05),
                c(0.08,0.05,0.4))
    pi2 <- c(0.5,0.5)
    pi3 <- c(0.4,0.4,0.2)
    vars2 <- vector('list',length(ns))
    vars3 <- vector('list',length(ns))
    names(vars2) <- ns
    names(vars3) <- ns
    for(i in 1:length(ns)){
        n <- ns[i]
        cat(n,"\n")
        ng2 <- n*pi2
        ng3 <- n*pi3
        cls2 <- rep(1:2,times=ng2)
        cls3 <- rep(1:3,times=ng3)
        vars2[[i]] <- vector('list',2)
        vars2[[i]][[1]] <- vector('list',N)
        vars2[[i]][[2]] <- vector('list',N)
        vars3[[i]] <- vector('list',3)
        vars3[[i]][[1]] <- vector('list',N)
        vars3[[i]][[2]] <- vector('list',N)
        vars3[[i]][[3]] <- vector('list',N)
        for(j in 1:N){
            set.seed(seed+i*N*2+j)
            g2 <- sample_sbm(n,B2,ng2,directed=FALSE)
            emb2 <- ase(g2,verbose=FALSE,max.d=20,scale.by.values=FALSE)
            vars2[[i]][[1]][[j]] <- var(emb2[which(cls2==1),])
            vars2[[i]][[2]][[j]] <- var(emb2[which(cls2==2),])
            g3 <- sample_sbm(n,B3,ng3,directed=FALSE)
            emb3 <- ase(g3,verbose=FALSE,max.d=20,scale.by.values=FALSE)
            vars3[[i]][[1]][[j]] <- var(emb3[which(cls3==1),])
            vars3[[i]][[2]][[j]] <- var(emb3[which(cls3==2),])
            vars3[[i]][[3]][[j]] <- var(emb3[which(cls3==3),])
        }
    }

    save(vars2,vars3,N,ns,seed,
         B2,B3,pi2,pi3,
         file='fig5.RData')
}

## Fig 5
figure5 <- function()
{
    files1 <- files %>% filter(grepl(glob2rx("fig5*"), name))
    file <- drive_download(files1, overwrite = TRUE, verbose = FALSE)
    load(file$local_path)
#    load('fig5.RData')

    data <- vector('list',2*length(ns))
    names(data) <- c(rbind(ns,rep("",length(ns))))
    k <- 1
    for(i in 1:length(ns)){
        data[[k]] <- unlist(lapply(vars2[[i]][[1]],function(y) {
            z <- y[-(1:2),-(1:2)]
            c(z[upper.tri(z)])
        }))
        data[[k+1]] <- unlist(lapply(vars2[[i]][[2]],function(y) {
            z <- y[-(1:2),-(1:2)]
            c(z[upper.tri(z)])
        }))
        k <- k+2
    }
    m <- par('mar')
    par(mar=c(m[1],m[1],m[2:3]))
    boxplot(data,col=c('white','gray'),xlab="n",ylab=expression(Sigma[ij]),
            cex.lab=1.5)
    abline(h=0,lty=2)
#    dev.print(device=pdf,file='offcov.pdf') # Fig 5(a)
    a <- boxplot(data,col=c('white','gray'),outline=FALSE,notch=TRUE)
    boxplot(data,col=c('white','gray'),outline=FALSE,notch=TRUE,
            ylim=c(min(a$stats[2,]),max(a$stats[4,])),
            xlab="n",ylab=expression(Sigma[ij]),cex.lab=1.5)
    abline(h=0,lty=2)
#    dev.print(device=pdf,file='offcov2.pdf') # Fig 5(b)
    par(mar=m)
}

experiment6 <- function(N=100,midp=0.2,maxp=0.5,overlap=.1,n=500,
                        ng=rep(1/2,2),seed=43523,
                        p1="exp4_all.pdf",
                        p2="exp42.pdf")
{
    m <- 2
    means <- matrix(0,nrow=N,ncol=4)
    offdiag <- rep(0,N)
    ondiag <- matrix(0,nrow=N,ncol=4)
    for(i in 1:N){
        set.seed(i+seed)
        P <- matrix(0,nrow=m,ncol=m)
        P[lower.tri(P)] <- runif(choose(m,2),0,midp)
        P[upper.tri(P)] <- t(P[lower.tri(P)])
        diag(P) <- runif(m,midp-overlap,maxp)
        NG <- n*ng
        while(sum(NG)>n) ng[which.max(ng)] <- ng[which.max(ng)] - 1
        while(sum(NG)<n) ng[which.min(ng)] <- ng[which.min(ng)] + 1
        g <- sample_sbm(n,P,NG,directed=FALSE)
        cls <- rep(1:m,times=NG)
        emb <- ase(g,max.d=m+2)
        means[i,] <- apply(emb,2,mean)
        v <- var(emb[,3:4])
        offdiag[i] <- v[1,2]
        ondiag[i,] <- diag(var(emb))
    }
    boxplot(cbind(means,offdiag,ondiag,ondiag[,3]-ondiag[,4]),
            notch=TRUE,
            #las=2,
            names=c(expression(mu[1]),
                    expression(mu[2]),
                    expression(mu[3]),
                    expression(mu[4]),
                    expression(Sigma[34]),
                    expression(Sigma[11]),
                    expression(Sigma[22]),
                    expression(Sigma[33]),
                    expression(Sigma[44]),
                    expression(Sigma[33]-Sigma[44])))
    abline(h=0,lty=2)
    dev.print(device=pdf,file=p1)
    par(mfrow=c(2,1))
    boxplot(means,
            notch=TRUE,
            #las=2,
            ylab="Means",
            names=c(expression(mu[1]),
                    expression(mu[2]),
                    expression(mu[3]),
                    expression(mu[4])))
    abline(h=0,lty=2)
    abline(v=2.5,lty=2)
    boxplot(cbind(offdiag,ondiag,ondiag[,3]-ondiag[,4]),
            ylab="Variances",
            names=c(expression(Sigma[34]),
                    expression(Sigma[11]),
                    expression(Sigma[22]),
                    expression(Sigma[33]),
                    expression(Sigma[44]),
                    expression(Sigma[33]-Sigma[44])))
    abline(h=0,lty=2)
    dev.print(device=pdf,file=p2)
    par(mfrow=c(1,1))
    list(means=means,off=offdiag,on=ondiag,midp=midp,maxp=maxp,overlap=overlap,n=n,
         ng=ng,seed=seed,p1=p1,p2=p2)
}

## Fig 6
figure6 <- function()
{
    out2 <- experiment6(N=10000,midp=0.1,maxp=0.5,overlap=-.1,n=500,
                        ng=rep(1/2,2),seed=4353,
                        p1="exp4_all_2.pdf",
                        p2="exp42_2.pdf") # Fig 6
}

experiment9 <- function(p=0.115,P,
                        n=500,
                        N=100,
                        seed=0,
                        pg=c(0.5,0.5),
                        adj.diag=FALSE)
{
    if(missing(P)){
        P <- rbind(c(0.2,p),c(p,0.1))
    }
    ng <- pg*n

    k <- nrow(P)
    outS <- vector('list',N)
    for(i in 1:N){
        set.seed(i+seed)
        g <- sample_sbm(n,P,ng,directed=FALSE)

        outS[[i]] <- aseCluster(g,verbose=FALSE,max.k=6,max.d=6,
                                scale.by.values=TRUE,short.circuit=TRUE,adjust.diag=adj.diag)
        emb <- ase(g,max.d=6,scale.by.values=TRUE)
        memb1 <- predict(outS[[i]]$best.model,emb)$classification
        memb2 <- predict(outS[[i]]$reduced.model,
                         emb[,1:outS[[i]]$d])$classification
        outS[[i]]$ari <- adjustedRandIndex(rep(1:k,times=ng),memb1)
        outS[[i]]$ariR <- adjustedRandIndex(rep(1:k,times=ng),memb2)
    }
    outS
}

experiment9.2 <- function(p=0.115,P,
                          n=500,
                          N=100,
                          seed=0,
                          pg=c(0.5,0.5),
                          adj.diag=FALSE)
{
    if(missing(P)){
        P <- rbind(c(0.2,p),c(p,0.1))
    }
    ng <- pg*n

    k <- nrow(P)
    outES <- vector('list',N)
    for(i in 1:N){
        set.seed(i+seed)
        g <- sample_sbm(n,P,ng,directed=FALSE)

        outES[[i]] <- try(aseClusterEM(g,verbose=FALSE,max.k=6,max.d=6,
                                       scale.by.values=TRUE,
                                       bail=TRUE,adjust.diag=adj.diag))
        if(!inherits(outES[[i]],'try-error')){
            emb <- ase(g,max.d=6,scale.by.values=TRUE)
            memb1 <- clusterMix(emb,outES[[i]]$best.model,d=6)
            memb2 <- clusterMix(emb[,1:outES[[i]]$best.model$d],
                                outES[[i]]$reduced.model,
                                d=outES[[i]]$best.model$d)
            outES[[i]]$ari <- adjustedRandIndex(rep(1:k,times=ng),memb1)
            outES[[i]]$ariR <- adjustedRandIndex(rep(1:k,times=ng),memb2)
        } else {
            outES[[i]] <- list(ari=NA,ariR=NA,bic=NA,BIC=NA,G=NA,d=NA)
        }
    }
    outES
}


# Fig 9
figure9 <- function(n=500,N=100,ps=seq(0.005,0.200,by=0.005))
{
#    if(!dir.exists(odir)) dir.create(odir,recursive=TRUE)
    ofile <- paste0("exp_N",N,"_n",n,".RData")

    # if (!file.exists(ofile)) {
    #     out <- mclapply(ps,function(p){
    #         experiment9(n=n,p=p,N=N)
    #     },mc.cores=30)
    #     outE <- mclapply(ps,function(p){
    #         experiment9.2(n=n,p=p,N=N)
    #     },mc.cores=30)
    #     save(out,outE,ps,n,file=ofile)
    # } else {
    #     load(ofile)
    # }

    files1 <- files %>% filter(grepl(glob2rx(ofile), name))
    file <- drive_download(files1, overwrite = TRUE, verbose = FALSE)
    load(file$local_path)


    a1 <- lapply(out,function(x) unlist(lapply(x,'[[','ari')))
    cols <- rep('white',length(ps))
    inds <- which(ps>=0.09 & ps<=0.115)
    cols[inds] <- 'gray'
    x1 <- boxplot(a1,names=ps,ylab="Adjusted Rand Index",xlab="p",
                  ylim=0:1,col=cols)
#    dev.print(device=pdf,file=paste0("Picts/aris",n,".pdf")) # Fig 9
}

figure1(K=2)
figure1(K=3)
figure2()
figure3()
figure5()
figure6()
figure9()

