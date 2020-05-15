library(parallel)
library(mvtnorm)
library(igraph)
library(RSpectra)
library(Matrix)
library(proxy)

ase <- function(g,verbose=FALSE,adjust.diag=FALSE,
                normalize=FALSE,scale.by.values=FALSE,
                vector='u',max.d=12)
{
   if(verbose) cat("Extracting adjacency matrix\n")
   A <- as_adjacency_matrix(g)
   if(adjust.diag){
      if(verbose) cat("Adjusting the diagonal of the adjacency matrix\n")
      n <- gorder(g)
      A <- A + Diagonal(x=degree(g)/(n-1))
   }
   if(normalize){
      if(verbose) cat("Normalizing the adjacency matrix\n")
      d <- degree(g)
      dinv <- sqrt(1/d)
      dinv[d==0] <- 0
      Dinv <- Diagonal(x=dinv)
      A <- Dinv %*% A %*% Dinv
   }
   if(verbose) cat("Computing the spectrum\n")
   if(!is.directed(g)){
      sv <- eigs(A,k=max.d)
      z <- sv$vectors
      if(scale.by.values) {
         ev <- sv$values
         ev[ev<=0] <- 1
         z <- scale(z,center=FALSE,scale=1/sqrt(ev))
      }
   } else {
      sv <- svds(A,k=max.d)
      if(tolower(vector)=='u'){
         z <- sv$u
         if(scale.by.values) z <- scale(z,center=FALSE,scale=1/sqrt(sv$d))
      } else if(tolower(vector)=='v'){
         z <- sv$v
         if(scale.by.values) z <- scale(z,center=FALSE,scale=1/sqrt(sv$d))
      } else {
         error("unknown vector specification")
      }
   }
   z
}

fig6 <- function(N=10000,out,
                 maxq=0.1,
				maxp=0.2,
				n=5000,
				d=20,
                ng=rep(1/2,2),
				outlines=1,
				seed=43523,
				savefile=paste0("RData/exp44_",N,"_",n,".RData"),
                pictfile=NULL)
{
	if(missing(out)){
		m <- 2
		means1 <- matrix(0,nrow=N,ncol=d)
		means2 <- matrix(0,nrow=N,ncol=d)
		ondiag1 <- matrix(0,nrow=N,ncol=d)
		ondiag2 <- matrix(0,nrow=N,ncol=d)
		bigoff1 <- vector('list',N)
		bigoff2 <- vector('list',N)
		simoff1 <- vector('list',N)
		simoff2 <- vector('list',N)
		dm12 <- rep(0,N)
		dmgt <- rep(0,N)
		off1_12 <- rep(0,N)
		off1_23 <- rep(0,N)
		off2_12 <- rep(0,N)
		off2_23 <- rep(0,N)
		for(i in 1:N){
		    cat("i = ", i, "\n")
			set.seed(i+seed)
			P <- matrix(0,nrow=m,ncol=m)
			P[lower.tri(P)] <- runif(choose(m,2),0,maxq)
			P[upper.tri(P)] <- t(P[lower.tri(P)])
			diag(P) <- sort(runif(m,max(P),maxp),decreasing=TRUE)
			NG <- n*ng
			g <- sample_sbm(n,P,NG,directed=FALSE)
			cls <- rep(1:m,times=NG)
			emb <- ase(g,max.d=d)
			means1[i,] <- apply(emb[cls==1,],2,mean)
			means2[i,] <- apply(emb[cls==2,],2,mean)
			dm12[i] <- proxy::dist(means1[i,1:2,drop=FALSE],means2[i,1:2,drop=FALSE])
			dmgt[i] <- proxy::dist(means1[i,-(1:2),drop=FALSE],means2[i,-(1:2),drop=FALSE])
			v1 <- var(emb[cls==1,])
			v2 <- var(emb[cls==2,])
			ondiag1[i,] <- diag(v1)
			ondiag2[i,] <- diag(v2)
			bigoff1[[i]] <- v1[upper.tri(v1[-(1:2),-(1:2)])]
			bigoff2[[i]] <- v2[upper.tri(v2[-(1:2),-(1:2)])]
			vvv <- diag(v1)[-(1:2)]
			s1 <- rmvnorm(N,mean=means1[i,-(1:2)],sigma=diag(vvv))
			sv1 <- var(s1)
			simoff1[[i]] <- sv1[upper.tri(sv1[-(1:2),-(1:2)])]

			vvv <- diag(v2)[-(1:2)]
			s2 <- rmvnorm(N,mean=means2[i,-(1:2)],sigma=diag(vvv))
			sv2 <- var(s2)
			simoff2[[i]] <- sv2[upper.tri(sv2[-(1:2),-(1:2)])]

			off1_12[i] <- v1[1,2]
			off2_12[i] <- v2[1,2]
			off1_23[i] <- v1[3,2]
			off2_23[i] <- v2[3,2]
		}
		bv1 <- unlist(bigoff1)
		bv2 <- unlist(bigoff2)
		sv1 <- unlist(simoff1)
		sv2 <- unlist(simoff2)
	   out <- list(m=m,N=N,
						maxq=maxq,
						maxp=maxp,
						n=n,
						d=d,
						ng=ng,
						seed=seed,
						means1=means1,
						means2=means2,
						ondiag1=ondiag1,
						ondiag2=ondiag2,
						bigoff1=bigoff1,
						bigoff2=bigoff2,
						simoff1=simoff1,
						simoff2=simoff2,
						dm12=dm12,
						dmgt=dmgt,
						off1_12=off1_12,
						off1_23=off1_23,
						off2_12=off2_12,
						off2_23=off2_23,
						bv1=bv1,bv2=bv2,
						sv1=sv1,sv2=sv2)
	   save(out,file=savefile)
	}
	with(out,{
		if(is.null(pictfile))
			pictfile <- paste0("figures/exp44_",N,"_",n,".pdf")
		layout(rbind(c(1,1,1,1,2,2,2,2),
						 c(1,1,1,1,2,2,2,2),
						 c(3,3,3,3,4,4,4,4),
						 c(3,3,3,3,4,4,4,4)))
		margin <- par('mar')
		tcl <- par('tcl')
		par(mar=c(4.1,4.1,3.1,1.5))

		boxplot(list(means1[,1],means1[,2],means2[,1],means2[,2]),
				  main="Informative Means",
				  notch=TRUE,
				  outline=1 %in% outlines,
				  names=c(expression(mu[1]^1),
							 expression(mu[2]^1),
							 expression(mu[1]^2),
							 expression(mu[2]^2)),
				  cex.names=1.5,
				  ylab="Mean",
				  xlab="Dimension")
		abline(h=0,lty=2)

		boxplot(list(c(means1[,-(1:2)]),c(means2[,-(1:2)])),
				  main="Redundant Means",
				  notch=TRUE,
				  outline=2 %in% outlines,
				  names=c(expression(mu[d>3]^1),expression(mu[d>3]^2)),
				  cex.names=1.5,
				  ylab="Mean",
				  xlab="Dimension")
		abline(h=0.0,lty=2)

		boxplot(list(ondiag1[,1],ondiag2[,1],
						 ondiag1[,2],ondiag2[,2]),
				  notch=TRUE,
				  outline=3 %in% outlines,
				  cex.names=1.5,
				  ylab="Informative Variance",
				  main="Variance",
				  names=c(
						  expression(sigma[1]^1),
						  expression(sigma[1]^2),
						  expression(sigma[2]^1),
						  expression(sigma[2]^2)))
		abline(h=0,lty=2)

		par(tcl=-0.025)
		boxplot(list(off1_12,off2_12,off1_23,off2_23,bv1,sv1,bv2,sv2),
				  notch=TRUE,
				  outline=4 %in% outlines,
				  main="Covariance",
				  cex.names=1.5,
				  col=c(rep('white',5),'gray','white','gray'),
				  ylab="Covariance",
				  names=c(
						  expression(Sigma[12]^1),
						  expression(Sigma[12]^2),
						  expression(Sigma[23]^1),
						  expression(Sigma[23]^2),
						  expression(Sigma[i!=j]^1),
						  expression(Sigma[i!=j]^1*s),
						  expression(Sigma[i!=j]^2),
						  expression(Sigma[i!=j]^2*s)))
		abline(h=0,lty=2)
		abline(v=2.5,lty=2)
		par(tcl=tcl)

		dev.print(device=pdf,file=pictfile)
		par(mfrow=c(1,1))
		par(mar=margin)
		return(invisible(out))
	})
}

