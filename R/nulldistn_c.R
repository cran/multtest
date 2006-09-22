#function to generate bootstrap null distribution
#theta0 is the value of the test statistics under the complete null hypthesis

boot.resample<-function(X,stat.closure,W=NULL,B=1000,theta0=0,tau0=1){	
	cat("running bootstrap...\n")
	Xnames<-dimnames(X)[[1]]
	X<-as.matrix(X)
	p<-nrow(X)
	n<-ncol(X)
	if(!(is.vector(W) | is.matrix(W) | is.null(W)))
		stop("W must be a vector or a matrix")	
	if(is.null(W))
		W<-matrix(1,nrow=p,ncol=n)
	if(is.vector(W)){
		if(length(W)==n)
			W<-matrix(W,nrow=p,ncol=n,byrow=TRUE)
		if(length(W)==p)
			W<-matrix(W,nrow=p,ncol=n)
		if(length(W)!=n & length(W)!=p)
			stop("Length of W does not match dim(X)")
	}
	if(is.matrix(W) & (dim(W)[1]!=p | dim(W)[2]!=n))
		stop("W and X must have same dimension")
	muboot<-matrix(0,nrow=p,ncol=B)
	samp<-sample(n,n*B,replace=TRUE)
	cat("iteration = ")
	muboot<-.C("bootloop",
           as.numeric(X),
           as.numeric(W),
           as.integer(p),
           as.integer(n),
           as.integer(B),
           muboot=double(p*B),
           as.integer(samp),
           body(stat.closure),
           environment(stat.closure),
	   NAOK=TRUE,	
           PACKAGE="multtest"
           )$muboot
	cat("\n")
	muboot<-matrix(muboot,nrow=p,ncol=B)
	dimnames(muboot)<-list(Xnames,paste(1:B))
	nas<-is.na(muboot)
	count<-0
	while(sum(nas)){
		count<-count+1
		if(count>1000)
			stop("Bootstrap null disrtibution computation 
terminating. Can not obtain distribution without missing values after 1000 
attempts.")
		nascols<-unique(col(muboot)[nas])
		for(b in nascols){
			samp<-sample(n,n,replace=TRUE)
			Xb<-X[,samp]
			Wb<-W[,samp]
			Tb<-get.Tn(Xb,stat.closure,Wb)
			muboot[,b]<-Tb[1,]/Tb[2,]
		}
		nas<-is.na(muboot)
	}
	(muboot-apply(muboot,1,mean)+theta0)*sqrt(pmin(1,tau0/apply(muboot,1,var)))
}



