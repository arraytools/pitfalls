# MC LI

rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

t.testv <- function(x, m1=ncol(x), m2=0, var.equal=TRUE) {
  # T-test by vectorization for two sample or one sample situation.
  # Missing values are denoted by NA and will be taken into account on the
  # degree of freedoms of the T-test.
  # depends on rowVars().
  #
  # Returns
  #  pval:
  #  t: class1 - class 2
  #  mean: class1 - class 2
  #  sd:
  #  n:
  #  df:
  m <- m1+m2
  if(!is.matrix(x)) x <- matrix(x, nr=1)
  if(ncol(x) != m) stop("m1 or m2 is not specified correctly")
  if(any(is.infinite(m))) stop("Some of expressions are +(-) Inf")
  pval <- rep(NA,nrow(x))
  if (m1<m) {
    # two sample t test
    n1 <- rowSums(!is.na(x[,1:m1,drop=FALSE]))
    n2 <- rowSums(!is.na(x[,(m1+1):m,drop=FALSE]))
    n <- n1+n2
    mean1 <- rowMeans(x[,1:m1,drop=FALSE],na.rm=TRUE)
    mean2 <- rowMeans(x[,(m1+1):m,drop=FALSE],na.rm=TRUE)
    if (var.equal) {
      # denom=pooled variance.
      denom <- (rowVars(x[,1:m1,drop=FALSE], na.rm=TRUE)*(n1-1) +
                  rowVars(x[,(m1+1):m,drop=FALSE], na.rm=TRUE)*(n2-1))/(n-2)
      tmp.sd <- sqrt(denom*((n1+n2)/(n1*n2)))
      tmp.df <- n-2
    }  else {
      # Behren's-Fisher problem
      denom <- rowVars(x[,1:m1,drop=FALSE], na.rm=TRUE)/n1 +
        rowVars(x[,(m1+1):m,drop=FALSE], na.rm=TRUE)/n2
      tmp.sd <- sqrt(denom)
      tmp.df <- denom^2 / ((rowVars(x[,1:m1,drop=FALSE], na.rm=TRUE))^2/(n1^2*(n1-1)) +
                             (rowVars(x[,(m1+1):m,drop=FALSE], na.rm=TRUE))^2/(n2^2*(n2-1)))
    }
    tpt <- (mean1-mean2)/tmp.sd
    ind <- which(n>2)
    pval[ind] <- 2 * (1 - pt(abs(tpt[ind]), tmp.df[ind]))
    result <- data.frame(pval=pval, t=tpt, mean=mean1-mean2, sd=tmp.sd, n=n, df=tmp.df)
  }
  else {
    # one sample t test
    n <- apply(x,1,function(x) sum(is.finite(x)))
    tmp.mean <- apply(x,1,mean,na.rm=TRUE)
    tmp.sd <- apply(x,1,sd,na.rm=TRUE)
    tmp <- which(n>1)
    tmp.t <- tmp.mean*sqrt(n)/tmp.sd # t-statistic
    pval[tmp] <- 2 * (1 - pt(abs(tmp.t[tmp]), n[tmp] - 1))
    result <- data.frame(pval=pval,t=tmp.t, mean=tmp.mean,sd=tmp.sd,n=n)
  }
  if (is.null(rownames(x)))
    rownames(result) <- paste("v",1:nrow(x), sep="")
  else
    rownames(result) <- rownames(x)
  result
}

ccp.train <- function(x, tt) {
  cc <- apply(x, 2, function(y) sum(tt * y))
  cc
}

ccp.predict <- function(c1, c2, tt, xnew) {
  cnew <- apply(xnew, 2, function(y) sum(tt * y))
  cmean <- (c1+c2)/2.0
  if (c1 <= c2) {
    pred <- ifelse(cnew <= cmean, 1, 2)
  } else {
    pred <- ifelse(cnew > cmean, 1, 2)
  }
  pred
}

rsbst <- function(p, pg) {
  n <- 20
  x <- matrix(rnorm(n*p), nr=p)
  out <- t.testv(x, 10, 10)
  indgene <- order(out$pval)[1:pg]
  ccp.tr <- ccp.train(x[indgene, ], out$t[indgene])
  ccp.pr <- ccp.predict(mean(ccp.tr[1:10]), mean(ccp.tr[11:20]),
                        out$t[indgene], x[indgene, ])
  err <- sum(abs(ccp.pr - c(rep(1, 10), rep(2, 10))))
  err
}

loocv1 <- function(p, pg) {
  # LOOCV after gene selection
  n <- 20
  x <- matrix(rnorm(n*p), nr=p)
  out <- t.testv(x, 10, 10)
  indgene <- order(out$pval)[1:pg]
  ccp.pr <- rep(NA, n)
  for(j in 1:n) {
    ccp.tr <- ccp.train(x[indgene, -j], out$t[indgene])
    if (j <= 10) {
      ccp.pr[j] <- ccp.predict(mean(ccp.tr[1:9]), mean(ccp.tr[10:19]),
                               out$t[indgene], x[indgene, j, drop = F])
    } else {
      ccp.pr[j] <- ccp.predict(mean(ccp.tr[1:10]), mean(ccp.tr[11:19]),
                               out$t[indgene], x[indgene, j, drop = F])
    }
  }
  err <- sum(abs(ccp.pr - c(rep(1, 10), rep(2, 10))))
  err
}

loocv2 <- function(p, pg) {
  # LOOCV before gene selection
  n <- 20
  x <- matrix(rnorm(n*p), nr=p)
  ccp.pr <- rep(NA, n)
  for(j in 1:n) {
    if (j <= 10) {
      n1 <- 9; n2 <- 10
    } else {
      n1 <- 10; n2 <- 9
    }
    out <- t.testv(x[, -j], n1, n2)
    indgene <- order(out$pval)[1:pg]
    ccp.tr <- ccp.train(x[indgene, -j], out$t[indgene])
    if (j <= 10) {
      ccp.pr[j] <- ccp.predict(mean(ccp.tr[1:9]), mean(ccp.tr[10:19]),
                               out$t[indgene], x[indgene, j, drop = F])
    } else {
      ccp.pr[j] <- ccp.predict(mean(ccp.tr[1:10]), mean(ccp.tr[11:19]),
                               out$t[indgene], x[indgene, j, drop = F])
    }
  }
  err <- sum(abs(ccp.pr - c(rep(1, 10), rep(2, 10))))
  err
}
