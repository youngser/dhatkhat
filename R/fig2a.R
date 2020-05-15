
figure1x <- function(max.d=80)
{
   load('fig1_2.RData')
   max.d <- min(max.d,ncol(mus2[[1]][[1]]))

	ylim <- c(-.0125,.0125)
   par(mfrow=c(2,2))
   boxplot(mus2[[1]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[1]),
	ylim=ylim,
   cex.lab=1.5)
   axis(2,at=c(-0.01,0,0.01))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=200")
   box()
   boxplot(mus2[[4]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[1]),
	ylim=ylim,
   cex.lab=1.5)
   axis(2,at=c(-0.01,0,0.01))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=2000")
   box()

   boxplot(mus2[[1]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[2]),
	ylim=ylim,
   cex.lab=1.5)
   axis(2,at=c(-0.01,0,0.01))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=200")
   box()
   boxplot(mus2[[4]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[2]),
	ylim=ylim,
   cex.lab=1.5)
   axis(2,at=c(-0.01,0,0.01))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=2000")
   box()
   dev.print(device=pdf,file='fig1meansx.pdf')

	ylim1 <- c(-.0125,.0125)
	ylim2 <- c(-.0003,.0003)
   par(mfrow=c(2,2))
   boxplot(mus2[[1]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[1]),
	ylim=ylim1,
   cex.lab=1.5)
   axis(2,at=c(-0.01,0,0.01))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=200")
   box()
   boxplot(mus2[[4]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[1]),
	ylim=ylim2,
   cex.lab=1.5)
   axis(2,at=c(-0.0002,0,0.0002))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=2000")
   box()

   boxplot(mus2[[1]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[2]),
	ylim=ylim1,
   cex.lab=1.5)
   axis(2,at=c(-0.01,0,0.01))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=200")
   box()
   boxplot(mus2[[4]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(mu[2]),
	ylim=ylim2,
   cex.lab=1.5)
   axis(2,at=c(-0.0002,0,0.0002))
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=2000")
   box()
   dev.print(device=pdf,file='fig1meansy.pdf')

   Y1 <- apply(rbind(mus2[[1]][[1]]^2,
                     mus2[[1]][[2]]^2)[,-(1:2)],2,mean)
   Y2 <- apply(rbind(mus2[[2]][[1]]^2,
                     mus2[[2]][[2]]^2)[,-(1:2)],2,mean)
   Y3 <- apply(rbind(mus2[[3]][[1]]^2,
                     mus2[[3]][[2]]^2)[,-(1:2)],2,mean)
   Y4 <- apply(rbind(mus2[[4]][[1]]^2,
                     mus2[[4]][[2]]^2)[,-(1:2)],2,mean)

   Y5 <- apply(rbind(mus2[[1]][[1]]^2,
                     mus2[[1]][[2]]^2)[,-(1:2)],2,median)
   Y6 <- apply(rbind(mus2[[2]][[1]]^2,
                     mus2[[2]][[2]]^2)[,-(1:2)],2,median)
   Y7 <- apply(rbind(mus2[[3]][[1]]^2,
                     mus2[[3]][[2]]^2)[,-(1:2)],2,median)
   Y8 <- apply(rbind(mus2[[4]][[1]]^2,
                     mus2[[4]][[2]]^2)[,-(1:2)],2,median)

   yl <- range(c(Y1,Y2,Y3,Y4))


   par(mfrow=c(1,1))

   plot(3:max.d,Y1,xlab="Dimension",ylab=expression(mu^2),ylim=yl,type='l')
   lines(Y2,col=2)
   lines(Y3,col=3)
   lines(Y4,col=4)
   b <- legend('topright',legend=paste0("n=",ns),col=1:4,lty=1)
   dev.print(device=pdf,file='fig1meanssq.pdf')
   lines(Y5,col=1,lty=2)
   lines(Y6,col=2,lty=2)
   lines(Y7,col=3,lty=2)
   lines(Y8,col=4,lty=2)
   legend(b$rect$left,b$rect$top-b$rect$h,c("Mean","Median"),lty=1:2)
   dev.print(device=pdf,file='fig1meanssqmed.pdf')

   par(mfrow=c(2,2))
   m <- par('mar')
   par(mar=c(m[1],m[2]+1,m[3:4]))
   boxplot(vars2[[1]][[1]][,-(1:2)],axes=FALSE,xlab='Dimension',ylab=expression(sigma[1]^2))
   axis(2)
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=200")
   box()
   boxplot(vars2[[4]][[1]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(sigma[1]^2))
   axis(2)
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=2000")
   box()

   boxplot(vars2[[1]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(sigma[2]^2))
   axis(2)
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=200")
   box()
   boxplot(vars2[[4]][[2]][,-(1:2)],names=3:max.d,axes=FALSE,xlab='Dimension',ylab=expression(sigma[2]^2))
   axis(2)
   axis(1,at=c(3,(1:(max.d/10))*10)-2,label=c(3,(1:(max.d/10))*10))
   title("n=2000")
   box()
   dev.print(device=pdf,file='fig1sigmas.pdf')
   par(mar=m)

   par(mfrow=c(2,2))
   boxplot(V[[1]][,-(1:2)],axes=FALSE,xlab='Dimension',ylab=expression(sigma^2))
   axis(2)
   axis(1,at=c(1,(1:(max.d/10))*10))
   title("n=200")
   boxplot(V[[2]][,-(1:2)],axes=FALSE,xlab='Dimension',ylab=expression(sigma^2))
   axis(2)
   axis(1,at=c(1,(1:(max.d/10))*10))
   title("n=500")
   boxplot(V[[3]][,-(1:2)],axes=FALSE,xlab='Dimension',ylab=expression(sigma^2))
   axis(2)
   axis(1,at=c(1,(1:(max.d/10))*10))
   title("n=1000")
   boxplot(V[[4]][,-(1:2)],axes=FALSE,xlab='Dimension',ylab=expression(sigma^2))
   axis(2)
   axis(1,at=c(1,(1:(max.d/10))*10))
   title("n=2000")

   d <- 7
   par(mfrow=c(2,2))
   m <- par('mar')
   par(mar=c(m[1],m[2]+1,m[3:4]))
   boxplot(vars2[[1]][[1]][,3:d],xlab='Dimension',ylab=expression(sigma[1]^2),
           names=3:d,notch=TRUE)
   title("n=200")
   box()
   boxplot(vars2[[4]][[1]][,3:d],xlab='Dimension',ylab=expression(sigma[1]^2),
           names=3:d,notch=TRUE)
   title("n=2000")
   box()

   boxplot(vars2[[2]][[2]][,3:d],xlab='Dimension',ylab=expression(sigma[2]^2),
           names=3:d,notch=TRUE)
   title("n=200")
   box()
   boxplot(vars2[[2]][[2]][,3:d],xlab='Dimension',ylab=expression(sigma[2]^2),
           names=3:d,notch=TRUE)
   title("n=2000")
   box()
   dev.print(device=pdf,file='fig1sigmasZ.pdf')
   par(mar=m)



   out1 <- vector('list',length(ns)*2)
   out2 <- vector('list',length(ns)*2)

   par(mfrow=c(2,2))
   k <- 1
   for(i in 1:length(ns)){
      out1[[k]] <- mus2[[i]][[1]][,3]
      out1[[k+1]] <- mus2[[i]][[2]][,3]
      out2[[k]] <- mus3[[i]][[1]][,3]
      out2[[k+1]] <- mus3[[i]][[2]][,3]
      k <- k+2
   }
   boxplot(out1,names=c(rbind(ns,rep("",length(ns)))),
           ylab=expression(mu),las=2,
           col=rep(c("white","lightgray"),times=length(ns)),
           notch=TRUE,main="d=3")

   k <- 1
   for(i in 1:length(ns)){
      out1[[k]] <- mus2[[i]][[1]][,4]
      out1[[k+1]] <- mus2[[i]][[2]][,4]
      out2[[k]] <- mus3[[i]][[1]][,4]
      out2[[k+1]] <- mus3[[i]][[2]][,4]
      k <- k+2
   }
   boxplot(out1,names=c(rbind(ns,rep("",length(ns)))),
           ylab=expression(mu),las=2,
           col=rep(c("white","lightgray"),times=length(ns)),
           notch=TRUE,main="d=4")
   k <- 1
   for(i in 1:length(ns)){
      out1[[k]] <- mus2[[i]][[1]][,5]
      out1[[k+1]] <- mus2[[i]][[2]][,5]
      out2[[k]] <- mus3[[i]][[1]][,5]
      out2[[k+1]] <- mus3[[i]][[2]][,5]
      k <- k+2
   }
   boxplot(out1,names=c(rbind(ns,rep("",length(ns)))),
           ylab=expression(mu),las=2,
           col=rep(c("white","lightgray"),times=length(ns)),
           notch=TRUE,main="d=5")

   k <- 1
   for(i in 1:length(ns)){
      out1[[k]] <- mus2[[i]][[1]][,6]
      out1[[k+1]] <- mus2[[i]][[2]][,6]
      out2[[k]] <- mus3[[i]][[1]][,6]
      out2[[k+1]] <- mus3[[i]][[2]][,6]
      k <- k+2
   }
   boxplot(out1,names=c(rbind(ns,rep("",length(ns)))),
           ylab=expression(mu),las=2,
           col=rep(c("white","lightgray"),times=length(ns)),
           notch=TRUE,main="d=6")

   par(mfrow=c(1,1))
   dev.print(device=pdf,file='mus3_6.pdf')

   x11(width=16,height=8)
   m <- par('mar')
   par(mfrow=c(2,1))
   par(mar=c(m[1],m[1],0.75,1))
   M1 <- lapply(mus2,function(x) (c(x[[1]][,-(1:2)])))
   M2 <- lapply(mus2,function(x) (c(x[[2]][,-(1:2)])))
   boxplot(M1,names=ns,notch=TRUE,xlab="n",ylab=expression(mu[1]),cex.lab=1.5)
   boxplot(M2,names=ns,notch=TRUE,xlab="n",ylab=expression(mu[2]),cex.lab=1.5)
   par(mar=m)
   dev.print(device=pdf,file='musall.pdf')

   par(mfrow=c(1,1))

}

