#main user-level function for empirical Bayes multiple hypothesis testing

EBMTP<-function(X,W=NULL,Y=NULL,Z=NULL,Z.incl=NULL,Z.test=NULL,na.rm=TRUE,test="t.twosamp.unequalvar",robust=FALSE,standardize=TRUE,alternative="two.sided",typeone="fwer",method="common.cutoff",k=0,q=0.1,alpha=0.05,smooth.null=FALSE,nulldist="boot.cs",B=1000,psi0=0,marg.null=NULL,marg.par=NULL,ncp=NULL,perm.mat=NULL,ic.quant.trans=FALSE,MVN.method="mvrnorm",penalty=1e-6,prior="conservative",bw="nrd",kernel="gaussian",seed=NULL,cluster=1,type=NULL,dispatch=NULL,keep.nulldist=TRUE,keep.rawdist=FALSE,keep.falsepos=FALSE,keep.truepos=FALSE,keep.errormat=FALSE,keep.Hsets=FALSE,keep.margpar=TRUE,keep.index=FALSE,keep.label=FALSE){
  ##sanity checks / formatting
  #X
  if(missing(X)) stop("Argument X is missing")
  if(inherits(X,"eSet")){ 
    if(is.character(Y)) Y<-pData(X)[,Y]
    if(is.character(Z)){
      if(Z%in%Y){
        Z<-Z[!(Z%in%Y)]
	warning(paste("Outcome Y=",Y,"should not be included in the covariates Z=",Z,". Removing Y from Z",sep=""))
	}
      Z<-pData(X)[,Z]
    }
    X<-exprs(X)
  }
  X<-as.matrix(X)
  dx<-dim(X)
  if(length(dx)==0) stop("dim(X) must have positive length")
  p<-dx[1]
  n<-dx[2]
  #W
  if(!is.null(W)){
    W[W<=0]<-NA
    if(is.vector(W) & length(W)==n) W <- matrix(rep(W,p),nr=p,nc=n,byrow=TRUE)
    if(is.vector(W) & length(W)==p) W <- matrix(rep(W,n),nr=p,nc=n)
    if(test%in%c("f","f.block","f.twoway","t.cor","z.cor")){
      warning("Weights can not be used with F-tests or tests of correlation parameters, arg W is being ignored.")
      W<-NULL
    }
  }
  #Y
  if(!is.null(Y)){
    if(is.Surv(Y)){
      if(test!="coxph.YvsXZ") stop(paste("Test ",test," does not work with a survival object Y",sep=""))
    }
    else{
      Y<-as.matrix(Y)
      if(ncol(Y)!=1) stop("Argument Y must be a vector")
    }
    if(nrow(Y)!=n) stop("Outcome Y has length ",nrow(Y),", not equal to n=",n)
  }
  #Z
  if(!is.null(Z)){
    Z<-as.matrix(Z)
    if(nrow(Z)!=n) stop("Covariates in Z have length ",nrow(Z),", not equal to n=",n,"\n")
    #Z.incl tells which columns of Z to include in model
    if(is.null(Z.incl)) Z.incl<-(1:ncol(Z))
    if(length(Z.incl)>ncol(Z)) stop("Number of columns in Z.incl ",length(Z.incl)," exceeds ncol(Z)=",ncol(Z))
    if(is.logical(Z.incl)) Z.incl<-(1:ncol(Z))[Z.incl]
    if(is.character(Z.incl) & length(Z.incl)!=sum(Z.incl%in%colnames(Z))) stop(paste("Z.incl=",Z.incl," names columns not in Z",sep=""))
    Za<-Z[,Z.incl]
    #Z.test tells which column of Z to test for an association
    if(test=="lm.XvsZ"){
      if(is.null(Z.test)){
        warning(paste("Z.test not specified, testing for association with variable in first column of Z:",colnames(Z)[1],sep=""))
	Z.test<-1
      }
      if(is.logical(Z.test)) Z.test<-(1:ncol(Z))[Z.test]
      if(is.character(Z.test) & !(Z.test%in%colnames(Z))) stop(paste("Z.test=",Z.test," names a column not in Z",sep=""))
      if(is.numeric(Z.test) & !(Z.test%in%(1:ncol(Z)))) stop("Value of Z.test must be >0 and <",ncol(Z))
      if(Z.test%in%Z.incl){
        Z.incl<-Z.incl[!(Z.incl%in%Z.test)]
	Za<-Z[,Z.incl]
      }
      Za<-cbind(Z[,Z.test],Za)
    }
    Z<-Za
    rm(Za)
  }
  #test
  TESTS<-c("t.onesamp","t.twosamp.equalvar","t.twosamp.unequalvar","t.pair","f","f.block","f.twoway","lm.XvsZ","lm.YvsXZ","coxph.YvsXZ","t.cor","z.cor")
  test<-TESTS[pmatch(test,TESTS)]
  if(is.na(test)) stop(paste("Invalid test, try one of ",TESTS,sep=""))

  #robust + see below with choice of nulldist
  if(test=="coxph.YvsXZ" & robust==TRUE)
    warning("No robust version of coxph.YvsXZ, proceeding with usual version")

  #alternative
  ALTS<-c("two.sided","less","greater")
  alternative<-ALTS[pmatch(alternative,ALTS)]
  if(is.na(alternative)) stop(paste("Invalid alternative, try one of ",ALTS,sep=""))

  #null values
  if(length(psi0)>1) stop(paste("In current implementation, all hypotheses must have the same null value. Number of null values: ",length(psi0),">1",sep=""))
  ERROR<-c("fwer","gfwer","tppfp","fdr")
  typeone<-ERROR[pmatch(typeone,ERROR)]
  if(is.na(typeone)) stop(paste("Invalid typeone, try one of ",ERROR,sep=""))
  if(any(alpha<0) | any(alpha>1)) stop("Nominal level alpha must be between 0 and 1.")
  nalpha<-length(alpha)
  reject<-
    if(nalpha) array(dim=c(p,nalpha),dimnames=list(rownames(X),paste("alpha=",alpha,sep="")))
    if(test=="z.cor" | test=="t.cor") matrix(nr=0,nc=0) # deprecated for correlations, rownames now represent p choose 2 edges - too weird and clunky in current state for output.
    else matrix(nr=0,nc=0)

  if(typeone=="fwer"){
    if(length(k)>1) k<-k[1]
    if(sum(k)!=0) stop("FWER control, by definition, requires k=0.  To control k false positives, please select typeone='gfwer'.")
  }

  if(typeone=="gfwer"){
    if(length(k)>1){
      k<-k[1]
      warning("Can only compute gfwer adjp for one value of k at a time (using first value). Use EBupdate() to get results for additional values of k.")
    }
    if(k==0) warning("gfwer(0) is the same as fwer.")
    if(k<0) stop("Number of false positives can not be negative.")
    if(k>=p) stop(paste("Number of false positives must be less than number of tests=",p,sep=""))
  }

  if(typeone=="tppfp"){
    if(length(q)>1){
      q<-q[1]
      warning("Can only compute tppfp adjp for one value of q at a time (using first value). Use EBupdate() to get results for additional values of q.")
    }
    if(q<0) stop("Proportion of false positives, q, can not be negative.")
    if(q>1) stop("Proportion of false positives, q, must be less than 1.")
  }

  #null distribution
  NULLS<-c("boot","boot.cs","boot.ctr","boot.qt","ic","perm")
  nulldist<-NULLS[pmatch(nulldist,NULLS)]
  if(is.na(nulldist)) stop(paste("Invalid nulldist, try one of ",NULLS,sep=""))
  if(nulldist=="perm") stop("EBMTP currently only available with bootstrap-based and influence curve null distribution methods.  One can, however, supply an externally created perm.mat for boot.qt marginal null distributions.")
  if(nulldist=="boot"){
    nulldist <- "boot.cs"
    warning("nulldist='boot' is deprecated and now corresponds to 'boot.cs'. Proceeding with default center and scaled null distribution.")
  }
  if(nulldist!="perm" & test=="f.block") stop("f.block test only available with permutation null distribution. Try test=f.twoway")
  if(nulldist=="ic" & keep.rawdist==TRUE) stop("Test statistics distribution estimation using keep.rawdist=TRUE is only available with a bootstrap-based null distribution")
  if(nulldist=="boot.qt" & robust==TRUE) stop("Quantile transform method requires parametric marginal nulldist.  Set robust=FALSE")
  if(nulldist=="boot.qt" & standardize==FALSE) stop("Quantile transform method requires standardized test statistics.  Set standardize=TRUE")
  if(nulldist=="ic" & robust==TRUE) stop("Influence curve null distributions available only for (parametric) t-statistics.  Set robust=FALSE")
  if(nulldist=="ic" & standardize==FALSE) stop("Influence curve null distributions available only for (standardized) t-statistics.  Set standardize=TRUE")
  if(nulldist=="ic" & (test=="f" | test=="f.twoway" | test=="f.block" | test=="coxph.YvsXZ")) stop("Influence curve null distributions available only for tests of mean, regression and correlation parameters. Cox PH also not yet implemented.")
  if(nulldist!="ic" & (test=="t.cor" | test=="z.cor")) stop("Tests of correlation parameters currently only implemented for influence curve null distributions")
  if((test!="t.cor" & test!="z.cor") & keep.index) warning("Matrix of indices only returned for tests of correlation parameters")

  ### specifically for sampling null test statistics with IC nulldist
  MVNS <- c("mvrnorm","Cholesky")
  MVN.method <- MVNS[pmatch(MVN.method,MVNS)]
  if(is.na(MVN.method)) stop("Invalid sampling method for IC-based MVN null test statistics.  Try either 'mvrnorm' or 'Cholesky'")

  #methods
  METHODS<-c("common.cutoff","common.quantile")
  method<-METHODS[pmatch(method,METHODS)]
  if(is.na(method)) stop(paste("Invalid method, try one of ",METHODS,sep=""))
  if(method=="common.quantile") stop("Common quantile procedure not currently implemented.  Common cutoff is pretty good, though.")

  #prior
  PRIORS<-c("conservative","ABH","EBLQV")
  prior<-PRIORS[pmatch(prior,PRIORS)]
  if(is.na(prior)) stop(paste("Invalid prior, try one of ",PRIORS,sep=""))

  #estimate 
  ftest<-FALSE
  if(test=="f" | test=="f.block"){
    ftest<-TRUE
    if(!is.null(W)) warning("Weighted F tests not yet implemented, proceding with unweighted version")
  }

    ##making a closure for the particular test
    theta0<-0
    tau0<-1
    stat.closure<-switch(test,
                         t.onesamp=meanX(psi0,na.rm,standardize,alternative,robust),
                         t.twosamp.equalvar=diffmeanX(Y,psi0,var.equal=TRUE,na.rm,standardize,alternative,robust),
                         t.twosamp.unequalvar=diffmeanX(Y,psi0,var.equal=FALSE,na.rm,standardize,alternative,robust),
                         t.pair={
                           uY<-sort(unique(Y))
                           if(length(uY)!=2) stop("Must have two class labels for this test")
                           if(trunc(n/2)!=n/2) stop("Must have an even number of samples for this test")
                           X<-X[,Y==uY[2]]-X[,Y==uY[1]]
                           meanX(psi0,na.rm,standardize,alternative,robust)
                         },
                         f={
                           theta0<-1
                           tau0<-2/(length(unique(Y))-1)
                           FX(Y,na.rm,robust)
                         },
                         f.twoway={
                           theta0<-1
                           tau0 <- 2/((length(unique(Y))*length(gregexpr('12', paste(Y, collapse=""))[[1]]))-1)
                           twowayFX(Y,na.rm,robust)
                         },
                         lm.XvsZ=lmX(Z,n,psi0,na.rm,standardize,alternative,robust),
                         lm.YvsXZ=lmY(Y,Z,n,psi0,na.rm,standardize,alternative,robust),
                         coxph.YvsXZ=coxY(Y,Z,psi0,na.rm,standardize,alternative),
                         t.cor=NULL,
                         z.cor=NULL)

    ##computing observed test statistics
    if(test=="t.cor" | test=="z.cor") obs<-corr.Tn(X,test=test,alternative=alternative,use="pairwise")
    else obs<-get.Tn(X,stat.closure,W)
    statistic <- (obs[3,]*obs[1,]/obs[2,]) #observed, with sign
    Tn <- obs[1,]/obs[2,]  # for sidedness, matching with mulldistn

      #Begin nulldists.  Permutation no longer included.
  if(nulldist=="boot.qt"){
  if(!is.null(marg.par)){
      if(is.matrix(marg.par)) marg.par <- marg.par
      if(is.vector(marg.par)) marg.par <- matrix(rep(marg.par,p),nr=p,nc=length(marg.par),byrow=TRUE)
        }
      if(is.null(ncp)) ncp = 0
      if(!is.null(perm.mat)){ 
        if(dim(X)[1]!=dim(perm.mat)[1]) stop("perm.mat must same number of rows as X.")
        }
    
      nstats <- c("t.twosamp.unequalvar","z.cor","lm.XvsZ","lm.YvsXZ","coxph.lmYvsXZ")
      tstats <- c("t.onesamp","t.twosamp.equalvar","t.pair","t.cor")
      fstats <- c("f","f.block","f.twoway")
      
      # If default, set values of marg.null to pass on. 
      if(is.null(marg.null)){
	  if(any(nstats == test)) marg.null="normal"
	  if(any(tstats == test)) marg.null="t"
	  if(any(fstats == test)) marg.null="f"
        }
      else{ # Check to see that user-supplied entries make sense.  
        MARGS <- c("normal","t","f","perm")
        marg.null <- MARGS[pmatch(marg.null,MARGS)]
        if(is.na(marg.null)) stop("Invalid marginal null distribution. Try one of: normal, t, f, or perm")
        if(any(tstats==test) & marg.null == "f") stop("Choice of test stat and marginal nulldist do not match")
        if(any(fstats==test) & (marg.null == "normal" | marg.null=="t")) stop("Choice of test stat and marginal nulldist do not match")
        if(marg.null=="perm" & is.null(perm.mat)) stop("Must supply a matrix of permutation test statistics if marg.null='perm'")
        if(marg.null=="f" & ncp < 0) stop("Cannot have negative noncentrality parameter with F distribution.")
      }
    
      # If default (=NULL), set values of marg.par. Return as m by 1 or 2 matrix.
      if(is.null(marg.par)){
		marg.par <- switch(test,
                          t.onesamp = n-1,
                          t.twosamp.equalvar = n-2,
                          t.twosamp.unequalvar = c(0,1),
                          t.pair = floor(n/2-1),
                          f = c(length(is.finite(unique(Y)))-1,dim(X)[2]- length(is.finite(unique(Y))) ),
                          f.twoway = {
                            c(length(is.finite(unique(Y)))-1, dim(X)[2]-(length(is.finite(unique(Y)))*length(gregexpr('12', paste(Y, collapse=""))[[1]]))-2)
                            },
                          lm.XvsZ = c(0,1),
                          lm.YvsXZ = c(0,1),
                          coxph.YvsXZ = c(0,1),
                          t.cor = n-2,
                          z.cor = c(0,1)
                          )
      marg.par <- matrix(rep(marg.par,dim(X)[1]),nr=dim(X)[1],nc=length(marg.par),byrow=TRUE)
              }
     else{ # Check that user-supplied values of marg.par make sense (marg.par != NULL)
       if((marg.null=="t" | marg.null=="f") & any(marg.par[,1]==0)) stop("Cannot have zero df with t or F distributions. Check marg.par settings")
       if(marg.null=="t" & dim(marg.par)[2]>1) stop("Too many parameters for t distribution.  marg.par should have length 1.")
       if((marg.null=="f" | marg.null=="normal") & dim(marg.par)[2]!=2) stop("Incorrect number of parameters defining marginal null distribution.  marg.par should have length 2.")
     }
}
  
    ##or computing influence curves
    if(nulldist=="ic"){
      rawdistn <- matrix(nr=0,nc=0)
      nulldistn<-switch(test,
                        t.onesamp=corr.null(X,W,Y,Z,test="t.onesamp",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat),
                        t.pair=corr.null(X,W,Y,Z,test="t.pair",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat),
                        t.twosamp.equalvar=corr.null(X,W,Y,Z,test="t.twosamp.equalvar",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat),
                        t.twosamp.unequalvar=corr.null(X,W,Y,Z,test="t.twosamp.unequalvar",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat),
                        lm.XvsZ=corr.null(X,W,Y,Z,test="lm.XvsZ",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat),
                        lm.YvsXZ=corr.null(X,W,Y,Z,test="lm.YvsXZ",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat),
                        t.cor=corr.null(X,W,Y,Z,test="t.cor",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat),
                        z.cor=corr.null(X,W,Y,Z,test="z.cor",alternative,use="pairwise",B,MVN.method,penalty,ic.quant.trans,marg.null,marg.par,perm.mat)
                        )
    }

    ## Cluster Checking
    if ((!is.numeric(cluster))&(!inherits(cluster,c("MPIcluster", "PVMcluster", "SOCKcluster"))))
       stop("Cluster argument must be integer or cluster object")
    ## Create cluster if cluster > 1 and load required packages on nodes
    if(is.numeric(cluster)){
      if(cluster>1){
    ## Check installation of packages
      have_snow <- qRequire("snow")
      if(!have_snow) stop("The package snow is required to use a cluster. Either snow is not installed or it is not in the standard library location.")
      if (is.null(type))
         stop("Must specify type argument to use a cluster. Alternatively, provide a cluster object as the argument to cluster.")
      if (type=="SOCK")
         stop("Create desired cluster and specify cluster object as the argument to cluster directly.")
      if ((type!="PVM")&(type!="MPI"))
         stop("Type must be MPI or PVM")
      else if (type=="MPI"){
         have_rmpi <- qRequire("Rmpi")
         if(!have_rmpi) stop("The package Rmpi is required for the specified type. Either Rmpi is not installed or it is not in the standard library location.")
      }
      else if (type=="PVM"){
         have_rpvm <- qRequire("rpvm")
         if(!have_rpvm) stop("The package rpvm is required for the specified type. Either rpvm is not installed or it is not in the standard library location.")
      }
      cluster <- makeCluster(cluster, type)
      clusterEvalQ(cluster, {library(Biobase); library(multtest)})
      if (is.null(dispatch)) dispatch=0.05
      }
    } else if(inherits(cluster,c("MPIcluster", "PVMcluster", "SOCKcluster"))){
      clusterEvalQ(cluster, {library(Biobase); library(multtest)})
      if (is.null(dispatch)) dispatch=0.05
    }

    ##computing the nonparametric bootstrap (null) distribution
    if(nulldist=="boot.cs" | nulldist=="boot.ctr" | nulldist=="boot.qt"){
      nulldistn<-boot.null(X,Y,stat.closure,W,B,test,nulldist,theta0,tau0,marg.null,marg.par,ncp,perm.mat,alternative,seed,cluster,dispatch,keep.nulldist,keep.rawdist)
     if(inherits(cluster,c("MPIcluster", "PVMcluster", "SOCKcluster")))  stopCluster(cluster)
    rawdistn <- nulldistn$rawboot
    nulldistn <- nulldistn$muboot
    }


    ##performing multiple testing
    #rawp values
    rawp<-apply((obs[1,]/obs[2,])<=nulldistn,1,mean)
    if(smooth.null & (min(rawp,na.rm=TRUE)==0)){
      zeros<-(rawp==0)
      if(sum(zeros)==1){
        den<-density(nulldistn[zeros,],to=max(obs[1,zeros]/obs[2,zeros],nulldist[zeros,],na.rm=TRUE),na.rm=TRUE)
	rawp[zeros]<-sum(den$y[den$x>=(obs[1,zeros]/obs[2,zeros])])/sum(den$y)
      }
      else{
        den<-apply(nulldistn[zeros,],1,density,to=max(obs[1,zeros]/obs[2,zeros],nulldistn[zeros,],na.rm=TRUE),na.rm=TRUE)
	newp<-NULL
	stats<-obs[1,zeros]/obs[2,zeros]
	for(i in 1:length(den)){
          newp[i]<-sum(den[[i]]$y[den[[i]]$x>=stats[i]])/sum(den[[i]]$y)
	}
        rawp[zeros]<-newp		
      }
      rawp[rawp<0]<-0
    }
    #c, cr, adjp - this is where the function gets a lot different from MTP.
    ### Begin nuts and bolts of EB here.

    ### Set G function of type I error rates
    error.closure <- switch(typeone, fwer=G.VS(V,S=NULL,tp=TRUE,bound=0),
                                     gfwer=G.VS(V,S=NULL,tp=TRUE,bound=k),
                                     tppfp=G.VS(V,S,tp=TRUE,bound=q),
                                     fdr=G.VS(V,S,tp=FALSE,bound=NULL)
                            )

    ### Generate guessed sets of true null hypotheses
    ### This function relates null and full densities.
    ### Sidedness should be accounted for above.
    H0.sets <- Hsets(Tn, nullmat=nulldistn, bw, kernel, prior, B, rawp=rawp) 
    EB.h0M <- H0.sets$EB.h0M
    prior.type <- prior
    prior.val <- H0.sets$prior
    lqv <- H0.sets$pn.out
    H0.sets <- H0.sets$Hsets.mat

    m <- length(Tn)

    ### B is defined in global environment.
    ### For adjusted p-values, just sort now and be able to get the index.
    ### We want to sort the test statistics in terms of their evidence against the null
    ### i.e., from largest to smallest.
    ord.Tn <- order(Tn,decreasing=TRUE)
    sort.Tn <- Tn[ord.Tn]
    Z.nulls <- nulldistn[ord.Tn,]*H0.sets[ord.Tn,]
    Tn.mat <- (1-H0.sets[ord.Tn,])*matrix(rep(sort.Tn,B),nr=m,nc=B)

    ### Rather than using a sieve of candidate cutoffs, for adjp, test statistics
    ### are used as cutoffs themselves.
    cutoffs <- sort.Tn
    clen <- m
    cat("counting guessed false positives...", "\n")
    Vn <- .Call(VScount,as.numeric(Z.nulls),as.numeric(cutoffs),as.integer(m),as.integer(B),as.integer(clen),NAOK=TRUE)
    cat("\n")
    Vn <- matrix(Vn, nr=clen, nc=B)

    if(typeone=="fwer" | typeone=="gfwer") Sn <- NULL
    else{
      cat("counting guessed true positives...", "\n")
      Sn <- .Call(VScount,as.numeric(Tn.mat),as.numeric(cutoffs),as.integer(m),as.integer(B),as.integer(clen),NOAK=TRUE)
      cat("\n")
      Sn <- matrix(Sn, nr=clen, nc=B)
    }

    G <-  error.closure(Vn,Sn)
    Gmeans <- rowSums(G,na.rm=TRUE)/B

    ### Now get adjps and rejection indicators.
    adjp <- rep(0,m)
    for(i in 1:m){
      adjp[i] <- min(Gmeans[i:m])
    }

    ### Now reverse order to go back to original order of test statistics.
    rev.order <- rep(0,m)
    for(i in 1:m){
      rev.order[i] <- which(sort.Tn==Tn[i])
    }
    adjp <- adjp[rev.order]
    if(keep.falsepos) Vn <- Vn[rev.order,]
    else Vn <- matrix(0,nr=0,nc=0)
    if(keep.truepos) Sn <- Sn[rev.order,]
    else Sn <- matrix(0,nr=0,nc=0)
    if(typeone=="fwer" | typeone=="gfwer") Sn <- matrix(0,nr=0,nc=0)
    if(keep.errormat) G <- G[rev.order,]
    else G <- matrix(0,nr=0,nc=0)
    if(!keep.Hsets) H0.sets <- matrix(0,nr=0,nc=0)
  
    # No confidence regions, but vector of rejections logicals, and cutoff, if applicable
    ### Generate matrix of rejection logicals.
    EB.reject <- reject
    if(test!="z.cor" & test!="t.cor") for(a in 1:nalpha) EB.reject[,a]<-adjp<=alpha[a]
    else EB.reject <- matrix(0,nr=0,nc=0)

    ### Grab test statistics corresponding to cutoff, based on adjp.
    #Leave out.
    #cutoff <- rep(0,nalpha)
    #for(a in 1:nalpha){
    #  if(sum(adjp<=alpha[a])>0){
    #   temp <- max(adjp[adjp<=alpha[a]])
    #   cutoff.ind <- which(adjp==temp)
    #   cutoff[a] <- max(Tn[cutoff.ind])
    # }
    #  else cutoff[a] <- NA
    #}

    #output results
    if(!keep.nulldist) nulldistn<-matrix(nr=0,nc=0)
    if(keep.rawdist==FALSE) rawdist<-matrix(nr=0,nc=0)
    if(is.null(Y)) Y<-matrix(nr=0,nc=0)
    if(nulldist!="boot.qt"){  
      marg.null <- vector("character")
      marg.par <- matrix(nr=0,nc=0)
    }
    if(!keep.label) label <- vector("numeric",0)
    if(!keep.index) index <- matrix(nr=0,nc=0)
    if(test!="z.cor" & test !="t.cor") index <- matrix(nr=0,nc=0)
    if(keep.index & (test!="z.cor" | test !="t.cor")){
      index <- t(combn(p,2))
      colnames(index) <- c("Var1","Var2")
    }
    names(adjp)<-names(rawp)
    out<-new("EBMTP",statistic=statistic,
      estimate=(if(ftest) vector("numeric",0) else obs[3,]*obs[1,]),
      sampsize=n,rawp=rawp,adjp=adjp,reject=EB.reject,
      rawdist=rawdistn,nulldist=nulldistn,nulldist.type=nulldist,
      marg.null=marg.null,marg.par=marg.par,
      label=label,falsepos=Vn,truepos=Sn,errormat=G,EB.h0M=EB.h0M,
      prior=prior.val,prior.type=prior.type,lqv=lqv,Hsets=H0.sets,
      index=index,call=match.call(),seed=as.integer(seed))

  return(out)
}



######################################################
######################################################
######################################################

### Function closure for different error rates.
G.VS <- function(V,S=NULL,tp=TRUE,bound){
  function(V,S){
    if(is.null(S)) g <- V     #FWER, GFWER
    else g <- V/(V+S)         #TPPFP, FDR
    if(tp==TRUE) {
      temp <- matrix(0,dim(g)[1],dim(g)[2])
      temp[g>bound] <- 1      #FWER, GFWER, TPPFP
      g <- temp
    }
    g
  }
}

### Adaptive BH estimate of the number of true null hypotheses.
ABH.h0 <- function(rawp){
  sortrawp <- sort(rawp)
  m <- length(rawp)
  ho.m <- rep(0,m)
  for(k in 1:length(rawp)){
    ho.m[k] <- (m+1-k)/(1-sortrawp[k])
  }
  grab <- min(which(diff(ho.m)>0))
  ho.hat <- ceiling(min(ho.m[grab],m))
  ho.hat
}

### Function for generating guessed sets via kernel density estimation.
### The marginal null is specified, although if boot.cs or boot.ctr is true,
### we can pool over the matrix of centered and scaled test statistics to
### estimate the null density.
### Also want the user to be able to set values of bw and kernel like they could
### using the density() function in R.
dens.est <- function(x,t,bw,kernel){
  dg <- density(t, from=x, to=x, bw=bw, kernel=kernel)
  dg$y[1]
}

Hsets <- function(Tn, nullmat, bw, kernel, prior, B, rawp){
### Full density estimation over vector of observed test statistics,
### saves on time and is asymptotic bootstrap distribution anyway.
### (As opposed to pooling over the whole matrix of raw tstats)
  f.Tn.est <- apply(as.matrix(Tn),1,dens.est,t=Tn, bw=bw, kernel=kernel)

  ### Obtain null density - use matrix of null test statistics...
  ### Ensures sidedness maintained more generally, especially in common-cutoff scenario
  dens.est.null <- approxfun(density(nullmat, bw=bw, kernel=kernel))
    f.Tn.null <- dens.est.null(Tn)

  ### pn represent local q-values obtained by density estimation
  ### Numbers might be so small and get returned NaN... true alts absolutely 
  pn <- pmin(1, f.Tn.null/f.Tn.est)
  pn[is.na(pn)] <- 0

  ### Do you want to relax the prior?
  if(prior=="conservative") priorval <- 1
  if(prior=="ABH")          priorval <- ABH.h0(rawp)/length(Tn)
  if(prior=="EBLQV")        priorval <- sum(pn,na.rm=TRUE)/length(Tn)

  pn.out <- pmin(1, priorval*pn)
  pn.out[is.na(pn.out)] <- 0

  # Draw Bernoullis for Ho matrix (Guessed sets of true null and true alternative hypotheses).
  # 1 = guessed true null, 0 = guessed true alternative hypotheses
  Hsets.mat <- matrix(rbinom(length(Tn)*B,1,pn.out),nr=length(Tn),nc=B)
  out <- list(Hsets.mat=Hsets.mat, EB.h0M=sum(pn,na.rm=TRUE)/length(Tn), prior=priorval, pn.out=pn.out)
  out
}

