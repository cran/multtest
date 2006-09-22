#functions that return a function (as in genefilter package)
# that takes a row of a data set as an arguement
# and returns a test statistic (possibly standardized)

#for one sample t
#paired t (where difference is the r.v.)
#Wilcoxon signed rank if robust=TRUE
meanX<-function(psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	function(x,w=NULL){
		if(is.null(w))
			w=rep(1,length(x))
		if(length(w)!=length(x))
			stop("x and w must have same length")
		x[!is.finite(w)]<-NA
		if(na.rm){
			drop<-is.na(x)|is.na(w)
			x<-x[!drop]
			w<-w[!drop]
		}
		if(robust) x<-(x>0)*rank(abs(x))
		n<-length(x)
		sumw<-sum(w)
		num<-sum(w*x)/sumw-psi0
		denom<-ifelse(is.na(num),NA,ifelse(standardize,sqrt(sum(w*(x-num)^2)/(sum(w*(1-(sum(w^2)/sumw^2))))),1))
		if(alternative=="two.sided"){
			snum<-sign(num)
		 	num<-abs(num)
		}
		if(alternative=="less"){
			snum<-(-1)
	 		num<-(-num)
		}
		if(alternative=="greater")
			snum<-1
		c(num,denom,snum)
	}
}

#two sample t (unequal variance (Welch) or equal variance)
#Wilcoxen rank sum if robust=TRUE and var.equal=FALSE
diffmeanX<-function(label,psi0=0,var.equal=FALSE,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	if(is.null(label))
		stop("A label variable is needed for this test")
	samp<-1:length(label)
	function(x,w=NULL){
		dep<-label[samp]
		if(is.null(w))
			w=rep(1,length(x))
		if(length(w)!=length(x))
			stop("x and w must have same length")
		x[!is.finite(w)]<-NA
		if(na.rm){
			drop<-is.na(x)|is.na(dep)|is.na(w)
			x<-x[!drop]
			xlabel<-dep[!drop]
			w<-w[!drop]
		}
		else
			xlabel<-dep
		if(robust) x<-rank(x)
		n<-length(x)
		lab1<-sort(unique(xlabel))[1]
		if(length(unique(xlabel))>2) 
			warning("More than 2 classes! Working with first unique label vs. rest.")
		sub1<-x[xlabel==lab1]
		sub2<-x[xlabel!=lab1]
		w1<-w[xlabel==lab1]
		w2<-w[xlabel!=lab1]
		m<-length(sub1)
		m1<-sum(w1*sub1)/sum(w1)
		m2<-sum(w2*sub2)/sum(w2)
		num<-m2-m1-psi0
		if(is.na(num))
			denom<-NA
		else{
			if(robust){
				df<-ifelse(sum(w)==1,sum(w),sum(w)-1)
				mm<-sum(w*x)/sum(w)
				denom<-ifelse(standardize,
					sqrt((1/m+1/(n-m))*sum(w*(x-mm)^2)/df),
					1)
			}
			else{
				if(standardize){
					df1<-sum(w1)-1
					df2<-sum(w2)-1
					df<-sum(w)-2
					denom<-ifelse(var.equal,
						sqrt((1/m+1/(n-m))*(sum(w1*(sub1-m1)^2)+sum(w2*(sub2-m2)^2))/df),
						sqrt((1/m*sum(w1*(sub1-m1)^2)/df1)+(1/(n-m)*sum(w2*(sub2-m2)^2)/df2)))
				}
				else 
					denom<-1
			}	
		}
		if(alternative=="two.sided"){
			snum<-sign(num)
		 	num<-abs(num)
		}
		if(alternative=="less"){
			snum<-(-1)
	 		num<-(-num)
		}
		if(alternative=="greater")
			snum<-1
		c(num,denom,snum)
	}
}

#F test for one way design
#Kruskal-Wallis rank.sum if robust=TRUE
FX<-function(label,na.rm=TRUE,robust=FALSE){
	if(is.null(label))
		stop("A label variable is needed for this test")
	samp<-1:length(label)
	function(x,w=NULL){
		dep<-label[samp]
		if(is.null(w)){
			if(na.rm){
				drop<-is.na(x)|is.na(dep)
				x<-x[!drop]
				xlabel<-dep[!drop]
			}
			else
				xlabel<-dep
			if(robust) x<-rank(x)
			n<-length(x)
			m<-mean(x)
			ulab<-sort(unique(xlabel))
			k<-length(ulab)
			nvec<-table(xlabel)
			mvec<-NULL
			for(i in 1:k){
				mvec[i]<-mean(x[xlabel==ulab[i]])
				x[xlabel==ulab[i]]<-x[xlabel==ulab[i]]-mvec[i]
			}
			num<-sum(nvec*(mvec-m)^2)
			denom<-sum(x^2)*(k-1)/(n-k)
		}
		else{
			if(length(w)!=length(x))
				stop("x and w must have same length")
			x[!is.finite(w)]<-NA
			if(na.rm){
				drop<-is.na(x)|is.na(dep)|is.na(w)
				x<-x[!drop]
				xlabel<-dep[!drop]
				w<-w[!drop]
			}
			else
				xlabel<-dep
			if(robust) x<-rank(x)
			n<-length(x)
#TODO: how to deal with weights in F?
			m<-mean(x)
			ulab<-sort(unique(xlabel))
			k<-length(ulab)
			nvec<-table(xlabel)
			mvec<-NULL
			for(i in 1:k){
				mvec[i]<-mean(x[xlabel==ulab[i]])
				x[xlabel==ulab[i]]<-x[xlabel==ulab[i]]-mvec[i]
			}
			num<-sum(nvec*(mvec-m)^2)
			denom<-sum(x^2)*(k-1)/(n-k)
		}
		c(num,denom,1)
	}
}

#could add standardize=FALSE and set denom=1 to return just the sum of squares.

#F statistic for block design with k treatments and l blocks
# The observations are ordered by block, and within 
#  each block, they are labeled using the integers 1 to k. 
#Friedman statistic if robust=TRUE
blockFX<-function(label,na.rm=TRUE,robust=FALSE){
	if(is.null(label))
		stop("A label variable is needed for this test")
	samp<-1:length(label)
	function(x,w=NULL){
		dep<-label[samp]
		if(is.null(w)){
			if(na.rm){
				drop<-is.na(x)|is.na(dep)
				x<-x[!drop]
				xlabel<-dep[!drop]
			}
			else
				xlabel<-dep
			n<-length(x)
			ulab<-sort(unique(xlabel))
			k<-length(ulab)
			l<-n/k
			if(round(l)*k!=n) stop("The blocks are not of equal size.")
			block<-sort(rep(1:l,k))[samp]
			ublock<-1:l
			if(robust){
			 	for(j in 1:l)
					x[block==j]<-rank(x[block==j])
			}
			m<-mean(x)			
			mlab<-mblock<-denom<-NULL
			for(i in 1:k){
				mlab[i]<-mean(x[xlabel==ulab[i]])
				for(j in 1:l){
					if(i==1) 
						mblock[j]<-mean(x[block==ublock[j]])
					denom[(i-1)*l+j]<-x[xlabel==ulab[i] & block==ublock[j]]-mlab[i]-mblock[j]+m
				}
			}
			num<-sum(l*(mlab-m)^2)
			denom<-sum(denom^2)/(l-1)
		}
		else{
			if(length(w)!=length(x))
				stop("x and w must have same length")
			x[!is.finite(w)]<-NA
			if(na.rm){
				drop<-is.na(x)|is.na(dep)|is.na(w)
				x<-x[!drop]
				xlabel<-dep[!drop]
				w<-w[!drop]
			}
			else
				xlabel<-dep
			n<-length(x)
			ulab<-sort(unique(xlabel))
			k<-length(ulab)
			l<-n/k
			if(round(l)*k!=n) stop("The blocks are not of equal size.")
			block<-sort(rep(1:l,k))[samp]
			ublock<-1:l
			if(robust){
			 	for(j in 1:l)
					x[block==j]<-rank(x[block==j])
			}
#TODO: how to deal with weights in block f?
			m<-mean(x)
			mlab<-mblock<-denom<-NULL
						for(i in 1:k){
				mlab[i]<-mean(x[xlabel==ulab[i]])
				for(j in 1:l){
					if(i==1)
						mblock[j]<-mean(x[block==ublock[j]])
					denom[(i-1)*l+j]<-x[xlabel==ulab[i] & block==ublock[j]]-mlab[i]-mblock[j]+m
				}
			}
			num<-sum(l*(mlab-m)^2)
			denom<-sum(denom^2)/(l-1)
		}
		c(num,denom,1)
	}
}

#could add standardize=FALSE and set denom=1 to return just the sum of squares.

#Z is a design *matrix*
# with variable of interest in first column
# and variables to adjust for in remaining columns
#Z is fixed for all columns of X
#gene expression is the outcome
lmX<-function(Z=NULL,n,psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	if(is.null(Z))
		Z<-matrix(1,n,1)
	else
		Z<-cbind(Z,rep(1,n))
	samp<-1:n
	function(x,w=NULL){
		covar<-Z[samp,]
		if(is.null(w)){
			if(na.rm){
				drop<-is.na(x)|apply(is.na(covar),1,sum)
				covar<-covar[!drop,]
				x<-x[!drop]
			}
			covar<-as.matrix(covar)
			if(robust){
				out<-rlm(covar,x)
				out$df.residual<-length(x)-out$rank
			}
			else
				out<-lm.fit(covar,x)
			denom<-ifelse(standardize,sqrt(sum(out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1]),1)
		}
		else{
			if(length(w)!=length(x))
				stop("x and w must have same length")
			x[!is.finite(w)]<-NA
			if(na.rm){
				drop<-is.na(x)|apply(is.na(covar),1,sum)|is.na(w)
				covar<-covar[!drop,]
				x<-x[!drop]
				w<-w[!drop]
			}
			covar<-as.matrix(covar)
			if(robust){
				out<-rlm(covar,x,w)
				out$df.residual<-length(x)-out$rank
			}	
			else
				out<-lm.wfit(covar,x,w)
			denom<-ifelse(standardize,sqrt(sum(w*out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1]),1)
		}
		num<-out$coef[1]-psi0
		if(alternative=="two.sided"){
			snum<-sign(num)
		 	num<-abs(num)
		}
		if(alternative=="less"){
			snum<-(-1)
	 		num<-(-num)
		}
		if(alternative=="greater")
			snum<-1
		c(num,denom,snum)
	}
}

#gene expression is the covariate
#y is an outcome of interest
#Z is any intercept or other covariates
#Z changes for each row of X
lmY<-function(Y,Z=NULL,n,psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	if(is.null(Y))
		stop("An outcome variable is needed for this test")
	if(length(Y)!=n)
		stop(paste("Dimension of outcome Y=",length(Y),", not equal dimension of data=",n,sep=""))
	if(is.null(Z))
		Z<-matrix(1,n,1)
	else
		Z<-cbind(Z,rep(1,n))
	samp<-1:n
	function(x,w=NULL){
		dep<-Y[samp]
		covar<-Z[samp,]
		covar<-cbind(x,covar)
		covar[!is.finite(w),]<-NA
		if(is.null(w)){
			if(na.rm){
				drop<-is.na(dep)|apply(is.na(covar),1,sum)
				covar<-covar[!drop,]
				xy<-dep[!drop]
			}
			else
				xy<-dep
			if(robust){
				out<-rlm(covar,xy)
				out$df.residual<-length(xy)-out$rank
			}
			else
				out<-lm.fit(covar,xy)
			denom<-ifelse(standardize,sqrt(sum(out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1]),1)
		}
		else{
			if(length(w)!=length(x))
				stop("x and w must have same length")
			if(na.rm){
				drop<-is.na(dep)|apply(is.na(covar),1,sum)|is.na(w)
				covar<-covar[!drop,]
				xy<-dep[!drop]
				w<-w[!drop]
			}
			else
				xy<-dep
			if(robust){
				out<-rlm(covar,xy,w)
				out$df.residual<-length(xy)-out$rank
			}
			else
				out<-lm.wfit(covar,xy,w)
			denom<-ifelse(standardize,sqrt(sum(w*out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1]),1)
		}
		num<-out$coef[1]-psi0
		if(alternative=="two.sided"){
			snum<-sign(num)
		 	num<-abs(num)
		}
		if(alternative=="less"){
			snum<-(-1)
	 		num<-(-num)
		}
		if(alternative=="greater")
			snum<-1
		c(num,denom,snum)
	}
}

#returns NA's if coxph fails
#strata is covariates to adjust for 
coxY<-function(surv.obj,strata=NULL,psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",init=NULL,method="efron"){
	autoload("coxph","survival")
    	if(!inherits(surv.obj,"Surv"))  #covers NULL case
        	stop("Response must be a survival object")
	if(!is.null(strata))
		strata<-as.matrix(strata)
	else
		strata<-strata(rep(1,nrow(surv.obj)))
    	samp<-1:nrow(surv.obj)
    	function(x,w=NULL){
        	if(!is.null(w)&length(w)!=length(x))
            		stop("x and w must have same length")
		dep<-surv.obj[samp,]
		covar<-strata[samp]
		if(na.rm){
			drop<-is.na(x)
			if(!is.null(w))
				drop<-drop+is.na(w)
			if(!is.null(strata))
				drop<-drop+is.na(covar)
			x<-x[!drop]
			w<-w[!drop]
			covar<-covar[!drop]
			dep<-dep[!drop,]
		}
		design<-cbind(x,rep(1,length(x)))		
		design[!is.finite(w),]<-NA
		control<-survival:::coxph.control()
		srvd<-try(survival:::coxph.fit(design,dep,strata=covar,init=init,control=control,weights=w,method=method,rownames=row.names(design)))
        	if(inherits(srvd,"try-error"))
	            return(c(NA,NA,NA))
	        denom<-ifelse(standardize,sqrt(srvd$var[1,1]),1)
       	 	num<-srvd$coef[1]-psi0
	        if(alternative=="two.sided"){
	            snum<-sign(num)
	            num<-abs(num)
	        }
	        if(alternative=="less"){  
	            snum<-(-1)
	            num<-(-num)       
	        }
	        if(alternative=="greater")
	            snum<-1
	        c(num,denom,snum)
	}   
}
	
#function that applies stat.closure to (X,W)
get.Tn<-function(X,stat.closure,W=NULL){
	wapply(X,1,stat.closure,W)	
}


