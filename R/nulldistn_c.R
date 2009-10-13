#functions to generate bootstrap null distribution
#theta0 is the value of the test statistics under the complete null hypthesis
#tau0 is the scaling parameter (upper  bound on variance of test statistics)

boot.null <- function(X,label,stat.closure,W=NULL,B=1000,test,nulldist,theta0=0,tau0=1,marg.null=NULL,marg.par=NULL,ncp=0,perm.mat,alternative="two.sided",seed=NULL,cluster=1,dispatch=0.05,keep.nulldist,keep.rawdist){
  cat("running bootstrap...\n")
  X<-as.matrix(X)
  n<-ncol(X)
  p<-nrow(X)
  if(!(is.vector(W) | is.matrix(W) | is.null(W))) stop("W must be a vector or a matrix")
  if(is.null(W))W<-matrix(1,nrow=p,ncol=n)
  if(is.vector(W)){
    if(length(W)==n) W<-matrix(W,nrow=p,ncol=n,byrow=TRUE)
    if(length(W)==p) W<-matrix(W,nrow=p,ncol=n)
    if(length(W)!=n & length(W)!=p) stop("Length of W does not match dim(X)")
  }
  if(is.matrix(W) & (dim(W)[1]!=p | dim(W)[2]!=n)) stop("W and X must have same dimension")

  # Dispatch to cluster
  if (is.numeric(cluster)) {
    if(!is.null(seed)) set.seed(seed)
    muboot <- boot.resample(X,label,p,n,stat.closure,W,B,test)
  }
  else {
    if(!is.null(seed)) clusterApply(cluster, seed, set.seed)
    else clusterApply(cluster, runif(length(cluster), max=10000000), set.seed)
    # Create vector of jobs to dispatch
          if ((dispatch > 0) & (dispatch < 1)){
            BtoNodes <- rep(B*dispatch, 1/dispatch)
          } else {
             BtoNodes <- rep(dispatch, B/dispatch)
          }
    FromCluster <- clusterApplyLB(cluster, BtoNodes, boot.resample,X=X,label=label,p=p,n=n,stat.closure=stat.closure,W=W, test=test)
    muboot <- matrix(unlist(FromCluster), nrow=nrow(X))
  }

  Xnames<-dimnames(X)[[1]]
  dimnames(muboot)<-list(Xnames,paste(1:B))

  #fill in any nas by resampling some more 
  nas<-(is.na(muboot)|muboot=="Inf"|muboot=="-Inf")
  count<-0
  while(sum(nas)){
    count<-count+1
    if(count>1000) stop("Bootstrap null distribution computation terminating. Cannot obtain distribution without missing values after 1000 attempts. This problem may be resolved if you try again with a different seed.")
    nascols<-unique(col(muboot)[nas])
    for(b in nascols){
      samp<-sample(n,n,replace=TRUE)
      Xb<-X[,samp]
      Wb<-W[,samp]
      if(p==1){
        Xb<-t(as.matrix(Xb))
	Wb<-t(as.matrix(Wb))
      }
      Tb<-get.Tn(Xb,stat.closure,Wb)
      muboot[,b]<-Tb[3,]*Tb[1,]/Tb[2,]
    }
    nas<-is.na(muboot)
  }

  rawboot <- matrix(nr=0,nc=0)
  if(keep.rawdist) rawboot <- muboot
  if(nulldist=="boot") muboot <- center.scale(muboot, theta0, tau0, alternative)
  if(nulldist=="boot.cs") muboot <- center.scale(muboot, theta0, tau0, alternative)
  if(nulldist=="boot.ctr") muboot <- center.only(muboot, theta0, alternative)
  if(nulldist=="boot.qt") muboot <- quant.trans(muboot, marg.null, marg.par, ncp, alternative, perm.mat)
  out <- list(muboot=muboot, rawboot=rawboot)
  out 

}

center.only <- function(muboot,theta0,alternative){
	muboot<-(muboot-apply(muboot,1,mean))+theta0
	if(alternative=="greater") muboot <- muboot
	else if(alternative=="less") muboot <- -muboot
	else muboot <- abs(muboot)
}

center.scale <- function(muboot, theta0, tau0, alternative){
  muboot<-(muboot-apply(muboot,1,mean))*sqrt(pmin(1,tau0/apply(muboot,1,var)))+theta0
  if(alternative=="greater") muboot <- muboot
  else if(alternative=="less") muboot <- -muboot
  else muboot <- abs(muboot)
}

quant.trans <- function(muboot, marg.null, marg.par, ncp, alternative, perm.mat){
### NB: Sanity checks occur outside this function at the beginning of MTP.
  m <- dim(muboot)[1]
  B <- dim(muboot)[2] 
  ranks <- t(apply(muboot,1,rank,ties.method="random"))
  Z.quant <- switch(marg.null,
                    normal = marg.samp(marg.null="normal",marg.par,m,B,ncp),
                    t = marg.samp(marg.null="t",marg.par,m,B,ncp),
                    f = marg.samp(marg.null="f",marg.par,m,B,ncp),
                    perm = perm.mat)
  Z.quant <- t(apply(Z.quant,1,sort))
### Left code like this for transparency. Could just as easily use quantile()
### for this first part, although it would be redundant.
  if(marg.null!="perm"){
    for(i in 1:m){
        Z.quant[i,] <- Z.quant[i,][ranks[i,]]
      }
    }
  else{
    Z.quant <- t(apply(Z.quant,1,quantile,probs=seq(0,1,length.out=B),na.rm=TRUE))
    for(i in 1:m){
        Z.quant[i,] <- Z.quant[i,][ranks[i,]]
      }
  }

  if(alternative=="greater") Z.quant <- Z.quant
  else if(alternative=="less") Z.quant <- -Z.quant
  else Z.quant <- abs(Z.quant)
  Z.quant
}


boot.resample <- function (X, label, p, n, stat.closure, W, B, test){
    muboot <- matrix(0, nrow = p, ncol = B)
    if (any(test == c("t.twosamp.equalvar", "t.twosamp.unequalvar",
        "f"))) {
        label <- as.vector(label)
        uniqlabs <- unique(label)
        num.group <- length(uniqlabs)
        groupIndex <- lapply(1:num.group, function(k) which(label ==
            uniqlabs[k]))
        obs <- sapply(1:num.group, function(x) length(groupIndex[[x]]))
        samp <- lapply(1:num.group, function(k) matrix(NA, nrow = B,
            ncol = obs[k]))
        for (j in 1:B) {
            for (i in 1:num.group) {
                uniq.obs <- 1
                count <- 0
                while (uniq.obs == 1) {
                  count <- count + 1
                  samp[[i]][j, ] <- sample(groupIndex[[i]], obs[i],
                    replace = TRUE)
                  uniq.obs <- length(unique(samp[[i]][j, ]))
                  if (count > 1000)
                    stop("Bootstrap null distribution computation terminating. Cannot obtain bootstrap sample with at least 2 unique observations after 1000 attempts. Sample size may be too small for bootstrap procedure but this problem may be resolved if you try again with a different seed.")
                }
            }
        }
        samp <- as.vector(t(matrix(unlist(samp), nrow = B, ncol = sum(obs))))
    }
    else if (test == c("f.twoway")) {
        label <- as.vector(label)
        utreat <- unique(label)
        num.treat <- length(utreat)
        num.block <- length(gregexpr("12", paste(label, collapse = ""))[[1]])
        ublock <- 1:num.block
        Breaks <- c(0, gregexpr(paste(c(num.treat, 1), collapse = ""),
            paste(label, collapse = ""))[[1]], n)
        BlockNum <- sapply(1:num.block, function(x) Breaks[x +
            1] - Breaks[x])
        block <- unlist(lapply(1:num.block, function(x) rep(x,
            BlockNum[x])))
        groupIndex <- lapply(1:num.block, function(j) sapply(1:num.treat,
            function(i) which(label == utreat[i] & block == ublock[j])))
        obs <- sapply(1:num.block, function(x) sapply(1:num.treat,
            function(y) length(groupIndex[[x]][[y]])))
        samp <- lapply(1:(num.treat * num.block), function(k) matrix(NA,
            nrow = B, ncol = obs[k]))
        for (k in 1:B) {
            for (i in 1:num.block) {
                for (j in 1:num.treat) {
                  uniq.obs <- 1
                  count <- 0
                  while (uniq.obs == 1) {
                    count <- count + 1
                    samp[[(i - 1) * num.treat + j]][k, ] <- sample(groupIndex[[i]][[j]],
                    obs[j, i], replace = TRUE)
                    uniq.obs <- length(unique(samp[[(i - 1) *
                      num.treat + j]][k, ]))
                    if (count > 1000)
                      stop("Bootstrap null distribution computation terminating. Cannot obtain bootstrap sample with at least 2 unique observations after 1000 attempts. Sample size may be too small for bootstrap procedure but this problem may be resolved if you try again with a different seed.")
                  }
                }
            }
        }
        samp <- as.vector(t(matrix(unlist(samp), nrow = B, ncol = sum(obs))))
      }
    else samp <- sample(n, n * B, replace = TRUE)
    cat("iteration = ")
    muboot <- .Call("bootloop", stat.closure, as.numeric(X),
        as.numeric(W), as.integer(p), as.integer(n), as.integer(B),
        as.integer(samp), NAOK = TRUE)
    cat("\n")
    muboot <- matrix(muboot, nrow = p, ncol = B)
}

