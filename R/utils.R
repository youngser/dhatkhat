
# calulate the BIC of the blockwise model 2, and then do clustering
# \mu_{d+1} = ... = \mu_{dmax} = 0; \sigma_{d+1} = ... = \sigma_{dmax}
# at least of init and k needs to be assigned
bicGmmBlock2 <- function(X, d, k, init = NULL) {
    # print("Block2")
    diff <- 1
    iter <- 0
    n <- dim(X)[1]
    dmax <- dim(X)[2]

    # initialization of parameters
    if (is.null(init)) {
        tao <- rep(1 / k, k)
        mu <- NULL
        Sigma1 <- NULL
        sigma2 <- NULL
        head <- 0
        for (i in 1 : k) {
            X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), 1 : d, drop = FALSE]
            X2 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), (d + 1) : dmax, drop = FALSE]
            head <- head + tao[i]
            mu[[i]] <- t(t(colMeans(X1)))
            Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
            for (l in 1 : dim(X1)[1]) {
                Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*%
                    (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
            }
            sigma2[[i]] <- 0
            for (l in 1 : dim(X2)[1]) {
                sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
            }
        }
    } else {
        label <- unique(init)
        numOfPoints <- rep(0, k)
        for (i in 1 : length(label)) {
            numOfPoints[i] <- sum(init == label[i])
        }
        indices <- sort.list(numOfPoints, decreasing = TRUE)
        if (length(label) >= k) {
            tao <- rep(0, k)
            mu <- NULL
            Sigma1 <- NULL
            sigma2 <- NULL
            for (i in 1 : k) {
                tao[i] <- sum(init == label[indices[i]])
                X1 <- X[init == label[indices[i]], 1 : d, drop = FALSE]
                X2 <- X[init == label[indices[i]], (d + 1) : dmax, drop = FALSE]
                mu[[i]] <- t(t(colMeans(X1)))
                Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
                for (l in 1 : dim(X1)[1]) {
                    Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*%
                        (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
                }
                sigma2[[i]] <- 0
                for (l in 1 : dim(X2)[1]) {
                    sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
                }
            }
            tao <- tao / sum(tao)
        } else {
            # length(label) < k
            tao <- rep(0, k)
            mu <- NULL
            Sigma1 <- NULL
            sigma2 <- NULL
            XX <- X[init == label[indices[1]], , drop = FALSE]
            numOfParti <- k - length(label) + 1
            numOfRows <- sum(init == label[indices[1]])
            for (i in 1 : numOfParti) {
                head <- floor(numOfRows * (i - 1) / numOfParti) + 1
                tail <- floor(numOfRows * (i) / numOfParti)
                tao[i] <- tail - head + 1
                X1 <- XX[head : tail, 1 : d, drop = FALSE]
                X2 <- XX[head : tail, (d + 1) : dmax, drop = FALSE]
                mu[[i]] <- t(t(colMeans(X1)))
                Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
                for (l in 1 : dim(X1)[1]) {
                    Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*%
                        (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
                }
                sigma2[[i]] <- 0
                for (l in 1 : dim(X2)[1]) {
                    sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
                }
            }
            for (i in (numOfParti + 1) : k) {
                tao[i] <- sum(init == label[indices[i - numOfParti + 1]])
                X1 <- X[init == label[indices[i - numOfParti + 1]], 1 : d, drop = FALSE]
                X2 <- X[init == label[indices[i - numOfParti + 1]], (d + 1) : dmax, drop = FALSE]
                mu[[i]] <- t(t(colMeans(X1)))
                Sigma1[[i]] <- matrix(rep(0, d*d), d, d)
                for (l in 1 : dim(X1)[1]) {
                    Sigma1[[i]] <- Sigma1[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*%
                        (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
                }
                sigma2[[i]] <- 0
                for (l in 1 : dim(X2)[1]) {
                    sigma2[[i]] <- sigma2[[i]] + sum(X2[l, ]^2) / (dim(X2)[1] * (dmax - d))
                }
            }
            tao <- tao / sum(tao)
        }
    }

    # Check if there is singularity
    for (j in 1 : k) {
        if (kappa(Sigma1[[j]]) > 1e10) {
            print("Singular at initialization")
            return (list(bic = 0,
                         classification = init,
                         d = d,
                         k = k,
                         iter = 0,
                         logliks = NULL,
                         tao = tao,
                         mu = mu,
                         Sigma1 = Sigma1,
                         sigma2 = sigma2))
        }
    }

    # initialization of f and t
    f <- matrix(rep(0, n * k), n, k)
    fTemp <- matrix(rep(0, n * k), n, k)
    t <- matrix(rep(0, n * k), n, k)
    for (i in 1 : n) {
        for (j in 1 : k) {
            f[i,j] <- tao[j] * (myDet(Sigma1[[j]])^(-0.5)) * sigma2[[j]]^(-(dmax - d) / 2) *
                exp(-(X[i, 1 : d, drop = FALSE] - t(mu[[j]])) %*%
                        solve(Sigma1[[j]]) %*% (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) / 2 -
                        sum(X[i, (d + 1) : dmax]^2) / (2 * sigma2[[j]]))
        }
    }
    loglik <- -(dmax * n / 2) * log(2 * pi)
    for (i in 1 : n) {
        fi <- sum(f[i,])
        if (fi > 1e-1000) {
            loglik <- loglik + log(fi)
        } else {
            loglik <- loglik - 1000
        }
    }
    #print(loglik)

    flag <- 0
    logliks <- rep(0, 100)
    logliks[iter + 1] <- loglik
    while ((flag == 0) && (diff > 1e-5) && (iter < 100)) {
        # print(iter)
        # if (iter == 8) {
        #   browser()
        # }
        # E step
        for (i in 1 : n) {
            if (sum(f[i,]) == 0) {
                for (j in 1 : k) {
                    t[i,j] <- 1/k
                }
            } else {
                for (j in 1 : k) {
                    t[i,j] <- f[i,j] / sum(f[i,])
                    if (t[i,j] < 1e-10) {
                        t[i,j] <- 1e-10
                    }
                }
            }
        }
        # M step
        for (j in 1 : k) {
            tao[j] <- sum(t[,j]) / n
            mu[[j]] <- matrix(rep(0, d), d, 1)
            for (i in 1 : n) {
                mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, 1 : d, drop = FALSE])
            }
            mu[[j]] <- mu[[j]] / sum(t[,j])
            Sigma1[[j]] <- matrix(rep(0, d^2), d, d)
            sigma2[[j]] <- 0
            for (i in 1 : n) {
                Sigma1[[j]] <- Sigma1[[j]] + t[i,j] * (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) %*%
                    (X[i, 1 : d, drop = FALSE] - t(mu[[j]]))
                sigma2[[j]] <- sigma2[[j]] + t[i,j] * sum(X[i, (d + 1) : dmax]^2)
            }
            Sigma1[[j]] <- Sigma1[[j]] / sum(t[,j])
            sigma2[[j]] <- sigma2[[j]] / ((dmax - d) * sum(t[,j]))
            # Calculate determinant and inverse of Sigma1 (181104)
        }
        # check convergence
        detSigma1 <- NULL
        invSigma1 <- NULL
        for (j in 1 : k) {
            if (kappa(Sigma1[[j]]) > 1e10) {
                print("Singular in the EM steps")
                numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2 + 1)
                #loglik <- preLoglik
                bic <- 2 * loglik - numOfParameters * log(n)
                classification <- rep(0, n)
                for (i in 1 : n) {
                    classification[i] <- which.max(f[i,])
                }
                return (list(bic = bic,
                             classification = classification,
                             d = d,
                             k = k,
                             iter = iter,
                             logliks = logliks,
                             tao = tao,
                             mu = mu,
                             Sigma1 = Sigma1,
                             sigma2 = sigma2))
            }
            detSigma1[[j]] <- myDet(Sigma1[[j]])
            invSigma1[[j]] <- solve(Sigma1[[j]])
        }

        for (i in 1 : n) {
            for (j in 1 : k) {
                if ((abs(detSigma1[[j]]) < 1e-1000) || (abs(sigma2[[j]]) < 1e-1000)) {
                    flag <- 1
                } else {
                    fTemp[i,j] <- tao[j] * ((detSigma1[[j]])^(-0.5)) * sigma2[[j]]^(-(dmax - d) / 2) *
                        exp(-(X[i, 1 : d, drop = FALSE] - t(mu[[j]])) %*%
                                invSigma1[[j]] %*% (t(X[i, 1 : d, drop = FALSE]) - mu[[j]]) / 2 -
                                sum(X[i, (d + 1) : dmax]^2) / (2 * sigma2[[j]]))
                }
            }
        }
        if (flag == 0) {
            f <- fTemp
        }
        preLoglik <- loglik
        loglik <- - (dmax * n / 2) * log(2 * pi)
        for (i in 1 : n) {
            fi <- sum(f[i,])
            if (fi > 1e-1000) {
                loglik <- loglik + log(fi)
            } else {
                loglik <- loglik - 1000
            }
        }
        #print(loglik)
        diff <- abs(loglik - preLoglik)
        iter <- iter + 1
        logliks[iter + 1] <- loglik
    }
    # calculate BIC
    numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2 + 1)
    bic <- 2 * loglik - numOfParameters * log(n)
    classification <- rep(0, n)
    for (i in 1 : n) {
        classification[i] <- which.max(f[i,])
    }
    return (list(bic = bic,
                 classification = classification,
                 d = d,
                 k = k,
                 iter = iter,
                 logliks = logliks,
                 tao = tao,
                 mu = mu,
                 Sigma1 = Sigma1,
                 sigma2 = sigma2))
}


# calculate the bic for VVV model with number of clusters k, then do clustering
bicGmmVVV <- function(X, k, tao = 0) {
    # print("VVV")
    diff <- 1
    iter <- 0
    n <- dim(X)[1]
    d <- dim(X)[2]

    # initialization of parameters
    if ((sum(tao) != 1) || (length(tao) != k)) {
        tao <- rep(1 / k, k)
    }
    mu <- NULL
    Sigma <- NULL
    head <- 0
    for (i in 1 : k) {
        X1 <- X[(1 + floor(n * head)) : (floor(n * (head + tao[i]))), , drop = FALSE]
        head <- head + tao[i]
        mu[[i]] <- t(t(colMeans(X1)))
        Sigma[[i]] <- matrix(rep(0, d*d), d, d)
        for (l in 1 : dim(X1)[1]) {
            Sigma[[i]] <- Sigma[[i]] + (t(X1[l, , drop = FALSE]) - mu[[i]]) %*%
                (X1[l, , drop = FALSE] - t(mu[[i]])) / dim(X1)[1]
        }
    }

    # Check if there is singularity
    for (j in 1 : k) {
        if (kappa(Sigma[[j]]) > 1e10) {
            print("Singular at initialization")
            return (list(bic = 0,
                         classification = init,
                         d = d,
                         k = k,
                         iter = 0,
                         logliks = NULL,
                         tao = tao,
                         mu = mu,
                         Sigma = Sigma))
        }
    }

    # initialization of f and t
    f <- matrix(rep(0, n * k), n, k)
    fTemp <- matrix(rep(0, n * k), n, k)
    t <- matrix(rep(0, n * k), n, k)
    for (i in 1 : n) {
        for (j in 1 : k) {
            if (abs(myDet(Sigma[[j]])) < 1e-10) {
                f[i,j] <- tao[j] * (1e-10)^(-0.5) *
                    exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*%
                            (1e10 * diag(d)) %*% (t(X[i, , drop = FALSE]) - mu[[j]]) / 2)
            } else {
                f[i,j] <- tao[j] * (myDet(Sigma[[j]])^(-0.5)) *
                    exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*%
                            solve(Sigma[[j]]) %*% (t(X[i, , drop = FALSE]) - mu[[j]]) / 2)
            }
        }
    }
    loglik <- -(d * n / 2) * log(2 * pi)
    for (i in 1 : n) {
        loglik <- loglik + log(sum(f[i,]))
    }

    flag <- 0
    while ((flag == 0) && (diff > 1e-5) && (iter < 100)) {
        # E step
        for (i in 1 : n) {
            if (sum(f[i,]) == 0) {
                for (j in 1 : k) {
                    t[i,j] <- 1/k
                }
            } else {
                for (j in 1 : k) {
                    t[i,j] <- f[i,j] / sum(f[i,])
                    if (t[i,j] < 1e-10) {
                        t[i,j] <- 1e-10
                    }
                }
            }
        }
        # M step
        for (j in 1 : k) {
            tao[j] <- sum(t[,j]) / n
            mu[[j]] <- matrix(rep(0, d), d, 1)
            for (i in 1 : n) {
                mu[[j]] <- mu[[j]] + t[i,j] * t(X[i, , drop = FALSE])
            }
            mu[[j]] <- mu[[j]] / sum(t[,j])
            Sigma[[j]] <- matrix(rep(0, d^2), d, d)
            for (i in 1 : n) {
                Sigma[[j]] <- Sigma[[j]] + t[i,j] * (t(X[i, , drop = FALSE]) - mu[[j]]) %*%
                    (X[i, , drop = FALSE] - t(mu[[j]]))
            }
            Sigma[[j]] <- Sigma[[j]] / sum(t[,j])
        }
        # check convergence
        for (i in 1 : n) {
            for (j in 1 : k) {
                if (abs(myDet(Sigma[[j]])) < 1e-10) {
                    flag <- 1
                } else {
                    fTemp[i,j] <- tao[j] * (myDet(Sigma[[j]])^(-0.5)) *
                        exp(-(X[i, , drop = FALSE] - t(mu[[j]])) %*%
                                solve(Sigma[[j]]) %*% (t(X[i, , drop = FALSE]) - mu[[j]]) / 2)
                }
            }
        }

        if (flag == 0) {
            f <- fTemp
        }

        preLoglik <- loglik
        loglik <- - (d * n / 2) * log(2 * pi)
        for (i in 1 : n) {
            loglik <- loglik + log(sum(f[i,]))
        }
        diff <- abs(loglik - preLoglik)
        iter <- iter + 1
    }
    # calculate BIC
    numOfParameters <- (k - 1) + k * (d + (d + 1) * d / 2)
    bic <- 2 * loglik - numOfParameters * log(n)
    classification <- rep(0, n)
    for (i in 1 : n) {
        classification[i] <- which.max(f[i,])
    }
    return (list(bic = bic,
                 loglik <- loglik,
                 classification = classification))
}

myDet <- function(X) {
    if (is.matrix(X)) {
        result <- det(X)
    } else {
        result <- X[1]
    }
    return (result)
}

getElbows <- function(dat, n = 3, threshold = FALSE, plot = F) {
    if (is.matrix(dat))
        d <- sort(apply(dat,2,sd), decreasing=T)
    else
        d <- sort(dat,decreasing=T)

    if (!is.logical(threshold))
        d <- d[d > threshold]

    p <- length(d)
    if (p == 0)
        stop(paste("d must have elements that are larger than the threshold ",
                   threshold), "!", sep="")

    lq <- rep(0.0, p)                     # log likelihood, function of q
    for (q in 1:p) {
        mu1 <- mean(d[1:q])
        mu2 <- mean(d[-(1:q)])              # = NaN when q = p
        sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) / (p - 1 - (q < p))
        lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
            sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
    }

    q <- which.max(lq)
    if (n > 1 && q < p) {
        q <- c(q, q + getElbows(d[(q+1):p], n=n-1, plot=FALSE))
    }

    if (plot==TRUE) {
        if (is.matrix(dat)) {
            sdv <- d # apply(dat,2,sd)
            plot(sdv,type="b",xlab="dim",ylab="stdev")
            points(q,sdv[q],col=2,pch=19)
        } else {
            plot(dat, type="b")
            points(q,dat[q],col=2,pch=19)
        }
    }

    return(q)
}

nonpsd.laplacian <- function(A)
{
    n <- nrow(A)
    s <- rowSums(A)/(n-1)
    diag(A) <- diag(A)+s
    return(A)
}

embed.graph <- function(A, dim, scaling = TRUE)
{
    suppressMessages(require(irlba))

    if (dim(A)[1]==dim(A)[2]) { ## diagonal augmentation
        L <- nonpsd.laplacian(A)
    } else {
        L <- A
    }

    if(is.null(dim)) {
        L.svd <- irlba(L)
        #        L.svd <- svd(L)
        #        dim <- scree.thresh(L.svd$d) # ignore this, use "getElbows" instead
    } else {
        L.svd <- irlba(L,dim,dim)
        #        L.svd <- svd(L)
    }

    L.svd.values <- L.svd$d[1:dim]
    L.svd.vectors <- L.svd$v[,1:dim]

    if(scaling == TRUE) { # projecting to sphere
        if(dim == 1)
            L.coords <- sqrt(L.svd.values) * L.svd.vectors
        else
            L.coords <- L.svd.vectors %*% diag(sqrt(L.svd.values))
    }
    else {
        L.coords <- L.svd.vectors
    }

    return(list(X=L.coords,D=L.svd.values))
}

plotAdjMat <- function(g)
{
    stopifnot(class(g) == "igraph")

    node_list <- get.data.frame(g, what = "vertices")
    if (ncol(node_list)==0) node_list$name <- 1:nrow(node_list)
    edge_list <- get.data.frame(g, what = "edges")
    all_nodes <- sort(node_list$name)

    plot_data <- edge_list %>% mutate(
        to = factor(to, levels = all_nodes),
        from = factor(from, levels = all_nodes))

    # Create the adjacency matrix plot
    ggplot(plot_data, aes(x = from, y = to)) +
        geom_raster() +
        theme_bw() +
        # Because we need the x and y axis to display every node,
        # not just the nodes that have connections to each other,
        # make sure that ggplot does not drop unused factor levels
        scale_x_discrete(drop = FALSE) +
        scale_y_discrete(drop = FALSE) +
        theme(
            # Rotate the x-axis lables so they are legible
            axis.text.x = element_text(angle = 90, hjust = 0),
            # Force the plot into a square aspect ratio
            aspect.ratio = 1,
            # Hide the legend (optional)
            legend.position = "none")
}


plotA <- function(g) {
    stopifnot(class(g)=="igraph")

    A <- as_adjacency_matrix(g)
    rownames(A) <- colnames(A) <- as.character(sprintf("%3d",1:nrow(A)))

    B <- list(
        from = rownames(A)[row(A)] %||% row(A),
        to = colnames(A)[col(A)] %||% col(A),
        value = A
    ) %>% map_dfc(as.vector)

    p <- ggplot(data=B, aes(x=from, y=to)) +
        coord_equal() +
        geom_tile(aes(fill=value), show.legend = FALSE) +
        scale_fill_gradientn(colors=gray(255:0/255), limit=c(0,1)) +
        theme(#axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              #strip.text.x = element_text(size=13),
              panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = "none") +
        labs(title="Adjaccency matrix", x="from", y="to")
    print(p)
}
## IRM stuff
p.R.z = function(R, z, prior.beta) {
    # Likelhood of data (R) given cluster assignment z. Evaluated by integrating
    # out possible values of the Bernoulli probabilities of the relations (which
    # is governed by prior.beta)
    z.uniq = unique(z)
    z.len = length(z.uniq)
    m.ab = matrix(0, nrow = z.len, ncol = z.len)
    m.ab.bar = matrix(0, nrow = z.len, ncol = z.len)
    # Loop through matrix. length(z.uniq))
    for (i in 1:nrow(R)) {
        for (j in 1:ncol(R)) {
            a = z[i]; b = z[j]
            R.ij = R[i, j]
            m.ab[a, b] = m.ab[a, b] + R.ij
            m.ab.bar[a, b] = m.ab.bar[a, b] + (1 - R.ij)
        }
    }

    ll = 0
    for (a in 1:nrow(m.ab)) {
        for (b in 1:nrow(m.ab)) {
            # ll <- foreach (a = 1:nrow(m.ab), .combine='c') %:%
            #   foreach (b = 1:nrow(m.ab), .combine='c') %dopar% {
            ll.ab = beta(m.ab[a, b] + prior.beta, m.ab.bar[a, b] + prior.beta) /
                beta(prior.beta, prior.beta)
            ll = ll + log(ll.ab)
            #  log(ll.ab)
        }
    }
    # sum(ll)
    ll
}

p.z = function(z, prior.gamma) {
    # Likelihood of this z being drawn from a CRP given GAMMA.
    # TODO: Sorting is probably really slow...
    z = sort(z)

    ll = 0
    # Keeping track of members already part of the clusters.
    # Requires that we know how many unique clusters there are
    n = rep(0, length(unique(z)))
    stopifnot(z[1] == 1)
    for (i in 1:length(z)) {
        a = z[i]
        n.a = n[a]
        if (n.a == 0) { # New cluster
            ll = ll + log(prior.gamma / (i - 1 + prior.gamma))
        } else { # Existing cluster
            ll = ll + log(n.a / (i - 1 + prior.gamma))
        }
        # Regardless, increment n.a
        n[a] = n[a] + 1
    }
    ll
}

renumber = function(arr) {
    # First, collapse extras
    arr = match(arr, unique(sort(arr)))
    # Second, how many clusters are there now?
    n.clus = length(unique(arr))
    a.map = rep(NA, n.clus)
    # Now loop from beinning to end of arr. When seeing a new number
    next.a = 1
    for (i in 1:length(arr)) {
        a = arr[i]
        if (is.na(a)) {
            next
        }
        if (is.na(a.map[a])) {
            a.map[a] = next.a
            arr[i] = next.a
            next.a = next.a + 1
        } else {
            arr[i] = a.map[a]
        }
    }
    arr
}

samp.z.probs = function(obs, prior.gamma) {
    # The probabilities of cluster assignments given already observed values.
    # The ith value is the length of obs + 1
    i = length(obs) + 1
    n.clus = length(unique(obs))
    n = rep(0, n.clus)
    for (j in 1:length(obs)) {
        a = obs[j]
        n[a] = n[a] + 1
    }
    # Now consider sampling i into any of the known as, OR a new a.
    p = rep(0, n.clus + 1)
    for (a in 1:n.clus) {
        # What is the probability that i belongs to a?
        n.a = n[a]
        stopifnot(n.a != 0)
        p[a] = (n.a / (i - 1 + prior.gamma))
    }
    # Sample from n.clus + 1 because you could sample a new cluster.
    p[n.clus + 1] = 1 - sum(p)
    p
}

samp.z = function(obs, prior.gamma) {
    # Samples from samp.z.probs.
    p = samp.z.probs(obs, prior.gamma)
    sample(length(p), 1, prob = p)
}

j = function(z.i, i, obs, prior.gamma) {
    # What was the loglik of sampling z.i as the ith position of obs?
    probs = samp.z.probs(renumber(obs[-i]), prior.gamma)
    log(probs[z.i])
}

# Simple 2-dimensional IRM. For each assignment "sample", I sweep through the
# cluster assignments of each individual point, and use a metropolis-hastings
# update.
irm = function(R, sweeps = 1000, prior.beta = 1, prior.gamma = 1, verbose = FALSE)
{

    Z = matrix(NA, nrow = sweeps, ncol = nrow(R))

    z = 1:nrow(R)
    n = length(z)

    # Sweeps.
    for (s in 1:sweeps) {
        for (i in 1:n) {
            if (verbose) cat("s = ", s, "/", sweeps, ", i = ", i, "/", n, "\n")
            # Propose changing this value by sampling from the existing CRP.
            # Specifically, since this is exchangeable, we imagine that z.i is the last
            # value to be sampled, and we sample from the probabilities.
            z.i = z[i]
            # Setup new vector. ith value to be filled in
            z.star = rep(z)
            z.star[i] = NA
            # Renumber (leaves NAs alone)
            z.star = renumber(z.star)

            # Get vector without the ith value and sample a new z.i
            z.i.star = samp.z(z.star[-i], prior.gamma)
            z.star[i] = z.i.star

            # Now compare likelihoods of 1) original z, and 2) new z.samp
            # Initial likelihood of original z and new z.samp
            log.ll = (p.R.z(R, z.star, prior.beta) +
                          p.z(z.star, prior.gamma)) -
                (p.R.z(R, z, prior.beta) + p.z(z, prior.gamma))
            # Next, relative likelihood of proposal distribution
            # likelihood of sampling z.i given z.i.star vs likelihood of sampling
            # z.i.star given z.i
            log.prop = j(z.i, i, z.star, prior.gamma) - j(z.i.star, i, z, prior.gamma)

            # Final acceptance ratio
            log.r = log.ll + log.prop

            if (log(runif(1)) < log.r) {
                z = renumber(z.star)
            }
        } # i
        # We ignore storing results of sweeps
        Z[s, ] = z
    } # s

    Z
}


# Collapse Z to strings
top.n = function(Z, n = 10) {
    Z.str = rep("", nrow(Z))
    for (i in 1:nrow(Z)) {
        Z.str[i] = paste(Z[i, ], collapse = "")
    }

    Z.counts = table(Z.str)
    # Only keep top 10 counts
    Z.top = head(sort(Z.counts, decreasing = TRUE), n = n)
    c(Z.top)
}

mode.irm = function(Z) {
    top.1 = names(top.n(Z, n = 1))
    as.numeric(strsplit(top.1, "")[[1]])
}

plot.R = function(R, assn = NULL) {
    if (is.null(assn)) {
        heatmap.2(R,
                  trace = 'none',
                  dendrogram = 'none',
                  Colv = 'none', Rowv = 'none',
                  key = FALSE,
                  col = c('#eeeeee', '#000000')
        )
    } else {
        assn.sort = sort(assn)
        # Where does assn change?
        assn.indices = c(1, 1 + which(diff(assn.sort) != 0)) - 1
        R.assn = R[order(assn), order(assn)]
        heatmap.2(R.assn,
                  trace = 'none',
                  dendrogram = 'none',
                  Colv = 'none', Rowv = 'none',
                  colsep = assn.indices,
                  rowsep = assn.indices,
                  key = FALSE,
                  col = c('#eeeeee', '#000000')
        )
    }
}
