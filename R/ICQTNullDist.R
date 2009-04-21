
# No robust correlation test statistics.
# Want to return a 3 by M matrix of observations.
corr.Tn <- function(X,test,alternative,use="pairwise"){
  P <- dim(X)[1]
  M <- P*(P-1)/2
  N <- dim(X)[2]
  VCM <- cov(t(X),use=use)
  Cor <- cov2cor(VCM)
  Cov.v <- VCM[lower.tri(VCM)] # vectorize.
  Cor.v <- Cor[lower.tri(Cor)] # vectorize.
  if(test=="t.cor") num <- sqrt(N-2)*Cor.v/sqrt(1-Cor.v^2)
  if(test=="z.cor") num <- sqrt(N-3)*0.5*log((1+Cor.v)/(1-Cor.v))
  denom <- 1
  if(alternative=="two.sided"){
			snum<-sign(num)
		 	num<-abs(num)
                      }
  else {
    if(alternative=="less"){
      snum<-(-1)
      num<-(-num)
    }
    else snum<-1
    }
  rbind(num,denom,snum)
}

ic.tests <- c("t.onesamp","t.pair","t.twosamp.equalvar","t.twosamp.unequalvar","lm.XvsZ","lm.YvsXZ","t.cor","z.cor")

corr.null <- function(X,W=NULL,Y=NULL,Z=NULL,test="t.twosamp.unequalvar",alternative="two-sided",use="pairwise",B=1000,MVN.method="mvrnorm",penalty=1e-6,ic.quant.trans=FALSE,marg.null=NULL,marg.par=NULL,perm.mat=NULL){
  # Most sanity checks conducted already...
  p <- dim(X)[1]
  m <- dim(X)[1] 
  n <- dim(X)[2] 
  cat("calculating vector influence curve...", "\n")

  if(test=="t.onesamp" | test=="t.pair"){
    #t.pair sanity checks and formatting done in stat.closure section
    #in test.R
    if(is.null(W)) IC.Cor <- cor(t(X),use=use)
    else IC.Cor <- IC.CorXW.NA(X,W,N=n,M=p,output="cor")
  }

  if(test=="t.twosamp.equalvar" | test=="t.twosamp.unequalvar"){
    uY<-sort(unique(Y))
    if(length(uY)!=2) stop("Must have two class labels for this test")
    n1 <- sum(Y==uY[1])
    n2 <- sum(Y==uY[2])
    if(is.null(W)){
      cov1 <- cov(t(X[,Y==uY[1]]),use=use)
      cov2 <- cov(t(X[,Y==uY[2]]),use=use)
    }
    else{
      cov1 <- IC.CorXW.NA(X[,Y==uY[1]],W[,Y==uY[1]],N=n1,M=p,output="cov")
      cov2 <- IC.CorXW.NA(X[,Y==uY[2]],W[,Y==uY[2]],N=n2,M=p,output="cov")
    }
    newcov <- cov1/n1 + cov2/n2
    IC.Cor <- cov2cor(newcov)
  }

  # Regression ICs written to automatically incorporate weights.
  # If W=NULL, then give equal weights.
  if(test=="lm.XvsZ"){
    if(is.null(Z)) Z <- matrix(1,nr=n,nc=1)
    else Z <- cbind(Z,1)
    if(is.null(W)) W <- matrix(1/n,nr=p,nc=n)
    IC.i <- matrix(0,nr=m,nc=n)
    for(i in 1:m){
      drop <- is.na(X[i,]) | is.na(rowSums(Z)) | is.na(W[i,])
      x <- as.numeric(X[i,!drop])
      z <- Z[!drop,]
      w <- W[i,!drop]
      nn <- n-sum(drop)
      EXtWXinv <- solve(t(z)%*%(w*diag(nn))%*%z)*sum(w)
      res.m <- lm.wfit(z,x,w)$res
      if(sum(drop)>0) res.m <- insert.NA(which(drop==TRUE),res.m)
      EXtWXinvXt <- rep(0,n)
      for(j in 1:n){
        EXtWXinvXt[j] <- (EXtWXinv%*%(t(Z)[,j]))[1]
      }
      IC.i[i,] <- res.m * EXtWXinvXt
    }
    IC.Cor <- IC.Cor.NA(IC.i,W,N=n,M=p,output="cor")
  }
  
  if(test=="lm.YvsXZ"){
    if(is.null(Y)) stop("An outcome variable is needed for this test")
    if(length(Y)!=n) stop(paste("Dimension of outcome Y=",length(Y),", not equal dimension of data=",n,sep=""))
    if(is.null(Z)) Z <- matrix(1,n,1)
    else Z <- cbind(Z,1)
    if(is.null(W)) W <- matrix(1,nr=p,nc=n)
    IC.i <- matrix(0,nr=m,nc=n)
    for(i in 1:m){
      drop <- is.na(X[i,]) | is.na(rowSums(Z)) | is.na(W[i,])
      x <- as.numeric(X[i,!drop])
      z <- Z[!drop,]
      w <- W[i,!drop]
      y <- Y[!drop]
      nn <- n-sum(drop)
      xz <- cbind(x,z)
      XZ <- cbind(X[i,],Z)
      EXtWXinv <- solve(t(xz)%*%(w*diag(nn))%*%xz)*sum(w)
      res.m <- lm.wfit(xz,y,w)$res
      if(sum(drop)>0) res.m <- insert.NA(which(drop==TRUE),res.m)
      EXtWXinvXt <- rep(0,n)
      for(j in 1:n){
        EXtWXinvXt[j] <- (EXtWXinv%*%(t(XZ)[,j]))[1]
      }
      IC.i[i,] <- res.m * EXtWXinvXt
    }
    IC.Cor <- IC.Cor.NA(IC.i,W,N=n,M=p,output="cor")
  }
  
  if(test=="t.cor" | test=="z.cor"){
    if(!is.null(W)) warning("Weights not currently implemented for tests of correlation parameters.  Proceeding with unweighted version")
    # Change of dimension
    P <- dim(X)[1] -> p # Number of variables.
    M <- P*(P-1)/2 -> m # Actual number of pairwise hypotheses.
    N <- dim(X)[2] -> m
    ind <- t(combn(P,2))
    VCM <- cov(t(X),use="pairwise")
    Cor <- cov2cor(VCM)
    Vars <- diag(VCM)
    Cov.v <- VCM[lower.tri(VCM)] # vectorize.
    Cor.v <- Cor[lower.tri(Cor)] # vectorize.
    X2 <- X*X
    EX <- rowMeans(X,na.rm=TRUE)
    E2X <- rowMeans(X2,na.rm=TRUE)
    Var1.v <- Vars[ind[,1]]
    Var2.v <- Vars[ind[,2]]
    EX1.v <- EX[ind[,1]]
    EX2.v <- EX[ind[,2]]
    E2X1.v <- E2X[ind[,1]]
    E2X2.v <- E2X[ind[,2]]
    X.vec1 <- X[ind[,1],]
    X.vec2 <- X[ind[,2],]
    X.vec12 <- X.vec1*X.vec2
    EX1X2.v <- rowMeans(X.vec12,na.rm=TRUE)

    cons <- 1/sqrt(Var1.v*Var2.v)
    gradient <- matrix(1,nr=M,nc=5)
    gradient[,1] <- EX1.v*Cov.v/Var1.v - EX2.v
    gradient[,2] <- EX2.v*Cov.v/Var2.v - EX1.v
    gradient[,3] <- -0.5*Cov.v/Var1.v
    gradient[,4] <- -0.5*Cov.v/Var2.v

    IC.i <- matrix(0, nr=M, nc=N)
    for(i in 1:N){
      diffs.i <- diffs.1.N(X[ind[,1],i], X[ind[,2],i], EX1.v, EX2.v, E2X1.v, E2X2.v, EX1X2.v)
      IC.M <- rep(0,M)
      for(j in 1:M){
        IC.M[j] <- gradient[j,]%*%diffs.i[,j]
      }
      IC.i[,i] <- IC.M
    }
    IC.i <- cons * IC.i
    IC.Cor <- IC.Cor.NA(IC.i,W=NULL,N=n,M=M,output="cor")
  }

  if(ic.quant.trans==FALSE) cat("sampling null test statistics...", "\n\n")
  else cat("sampling null test statistics...", "\n")
  
  if(MVN.method=="mvrnorm") nulldist <- t(mvrnorm(n=B,mu=rep(0,dim(IC.Cor)[1]),Sigma=IC.Cor))
  if(MVN.method=="Cholesky"){
    IC.chol <- t(chol(IC.Cor+penalty*diag(dim(IC.Cor)[1])))
    norms <- matrix(rnorm(B*dim(IC.Cor)[1]),nr=dim(IC.Cor)[1],nc=B)
    nulldist <- IC.chol%*%norms
  }
  if(ic.quant.trans==TRUE){
    cat("applying quantile transform...", "\n\n")
    if(is.null(marg.null)){
	marg.null <- "t"
     	if(test=="t.cor" | test=="z.cor" | test=="t.twosamp.equalvar") marg.par <- matrix(rep(dim(X)[2]-2,dim(IC.Cor)[1]),nr=dim(IC.Cor)[1],nc=1)
        if(test=="lm.XvsZ") marg.par <- matrix(rep(dim(X)[2]-dim(Z)[2],dim(IC.Cor)[1]),nr=dim(IC.Cor)[1],nc=1)
        if(test=="lm.YvsXZ")  marg.par <- matrix(rep(dim(X)[2]-dim(Z)[2]-1,dim(IC.Cor)[1]),nr=dim(IC.Cor)[1],nc=1)
        else marg.par <- matrix(rep(dim(X)[2]-1,dim(IC.Cor)[1]),nr=dim(IC.Cor)[1],nc=1)
    }
    if(test=="z.cor" & marg.null=="t") warning("IC nulldist for z.cor already MVN. Transforming to N-2 df t marginal distribution not advised.")
    if(marg.null!="t" & marg.null!="perm") stop("IC nulldists can only be quantile transformed to a marginal t-distribution or user-supplied marginal permutation distribution")
    if(marg.null=="t") nulldist <- tQuantTrans(nulldist,marg.null="t",marg.par,ncp=0,perm.mat=NULL)
    if(marg.null=="perm") nulldist <- tQuantTrans(nulldist,marg.null="perm",marg.par=NULL,ncp=NULL,perm.mat=perm.mat)
  }
  if(alternative=="greater") nulldist <- nulldist
  else if(alternative=="less") nulldist <- -nulldist
  else nulldist <- abs(nulldist)
  nulldist
}

# Function, given ICs for each individual, returns variance covariance
# matrix or corresponding correlation matrix.
IC.Cor.NA <- function(IC,W,N,M,output){
  n <- dim(IC)[2]
  m <- dim(IC)[1]
  if(is.null(W)){
    W <- matrix(1,nr=dim(IC)[1],nc=dim(IC)[2])
    Wnew <- W/rowSums(W,na.rm=TRUE) # Equal weight, NA handling.
  }
  else Wnew <- W/rowSums(W,na.rm=TRUE)
  IC.VC <- matrix(0,nr=m,nc=m)
  for(i in 1:n){
    temp <- crossprod(t(sqrt(Wnew[,i])*IC[,i]))
    temp[is.na(temp)] <- 0
    IC.VC <- IC.VC + temp
  }
  if(output=="cov") out <- IC.VC
  if(output=="cor") out <- cov2cor(IC.VC)
  out
}
 
# Weighted correlation. Generalizes cov.wt() to account for a matrix
# of weights. Uses IC formulation instead of sweep() and crossprod().
# May be slower/clunkier, but pretty transparent, and allows for NA
# handling much like cor(...,use="pairwise") would.  That is, each
# element of the correlation matrix returned uses the maximum amount
# of information possible in obtaining individual elements of that
# matrix.
IC.CorXW.NA <- function(X,W,N,M,output){
  n <- dim(X)[2]
  m <- dim(X)[1]
  XW <- X*W
  EXW <- rowSums(XW)/rowSums(W)
  ICW.i <- X-EXW
  Wnew <- W/rowSums(W,na.rm=T)
  IC.VC <- matrix(0,nr=m,nc=m)
  for(i in 1:n){
    temp <- crossprod(t(sqrt(Wnew[,i])*X[,i]))
    temp[is.na(temp)] <- 0
    IC.VC <- IC.VC + temp
  }
  if(output=="cov") out <- IC.VC
  if(output=="cor") out <- cov2cor(IC.VC)
  out
}

# For regression ICs, a function to insert NAs into appropriate locations
# of a vector of returned residuals.
insert.NA <- function(orig.NA, res.vec){
  for(i in 1:length(orig.NA)){
    res.vec <- append(res.vec, NA, after=orig.NA[i]-1)
  }
  res.vec
}

# For correlation ICS, a function to get diff vectors for all M.
# This is the difference between estimates for
# a sample size of one and a sample of size n.
diffs.1.N <- function(vec1, vec2, e1, e2, e21, e22, e12){
  diff.mat.1.N <- matrix(0,nr=5,nc=length(vec1))
  diff.mat.1.N[1,] <- vec1 - e1
  diff.mat.1.N[2,] <- vec2 - e2
  diff.mat.1.N[3,] <- vec1*vec1 - e21
  diff.mat.1.N[4,] <- vec2*vec2 - e22
  diff.mat.1.N[5,] <- vec1*vec2 - e12
  diff.mat.1.N
}

### For quantile transform, take a sample from the marginal null distribution.
marg.samp <- function(marg.null,marg.par,m,B,ncp){
out <- matrix(0,m,B)
for(i in 1:m){
  if(marg.null=="normal") out[i,] <- rnorm(B,mean=marg.par[i,1],sd=marg.par[i,2])
  if(marg.null=="t") out[i,] <- rt(B,df=marg.par[i,1],ncp)
  if(marg.null=="f") out[i,] <- rf(B,df1=marg.par[i,1],df2=marg.par[i,2],ncp)
}
out
}

### Quantile transform streamlined for IC nulldists.
tQuantTrans <- function(rawboot, marg.null, marg.par, ncp, perm.mat=NULL){
  m <- dim(rawboot)[1]
  B <- dim(rawboot)[2] 
  ranks <- t(apply(rawboot,1,rank,ties.method="random"))
  if(marg.null=="t") Z.quant <- marg.samp(marg.null="t",marg.par,m,B,ncp)
  if(marg.null=="perm") Z.quant <- perm.mat
  Z.quant <- t(apply(Z.quant,1,sort))
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
  Z.quant
}

### Effective df for two sample test of means, unequal var.
t.effective.df <- function(X,Y){
  uY<-sort(unique(Y))
  X1 <- X[Y==uY[1]]
  X2 <- X[Y==uY[2]]
  mu <- var(X2)/var(X1)
  n1 <- length(Y[Y==uY[1]])
  n2 <- length(Y[Y==uY[2]])
  df <- (((1/n1)+(mu/n2))^2)/(1/((n1^2)*(n1-1)) + (mu^2)/((n2^2)*(n2-1)))
  df
}






    












  
