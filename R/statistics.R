#for one sample t
#paired t (where difference is the r.v.)
#Wilcoxon signed rank if robust=TRUE
meanX<-function(psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	function(x,w=NULL, samp){
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
    if(is.na(num)) denom <- NA
    else{
     if(standardize) {denom <- sqrt(sum(w*(x-num)^2)/(sum(w*(1-(sum(w^2)/sumw^2)))))}
    else
       denom <- 1
     }
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
		c(num*sqrt(sumw),denom,snum)
	}
}

diffmeanX<-function(label,psi0=0,var.equal=FALSE,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	if(is.null(label))
		stop("A label variable is needed for this test")
	Samp<-1:length(label)
	function(x,w=NULL, samp=Samp){
		dep<-label[samp]
		if(is.null(w))
			w=rep(1,length(x))
		if(length(w)!=length(x))
			stop("x and w must have same length")
		x[!is.finite(w)]<-NA
		if(na.rm){
			drop<-is.na(x)|is.na(dep)|is.na(w)
			x<-x[!drop]
			xlabel<-as.vector(dep[!drop])
			w<-w[!drop]
		}

		else
			xlabel<-as.vector(dep)
    # Convert to 0,1 coding
    uniq <- sort.int(unique.default(xlabel))
 		if(length(uniq)>2)
			warning("More than 2 classes! Working with first unique label vs. rest.")
    lab1<-uniq[1]
    New0 <- which(xlabel==lab1)
    New1 <- which(xlabel!=lab1)
    xlabel <- as.numeric(replace(replace(xlabel, New1, 1), New0, 0))
    xlabel <- as.numeric(xlabel)

    # Check for at least 2 unique values in each group
    if(standardize & length(unique.default(x[xlabel==lab1]))==1) stop("Only one unique value in bootstrap sample for first group. Cannot calculate variance. This problem may be resolved if you try again with a different seed.")
    if(standardize & length(unique.default(x[xlabel!=lab1]))==1) stop("Only one unique value in bootstrap sample for second group. Cannot calculate variance. This problem may be resolved if you try again with a different seed.")
    n<-length(x)
    if(robust) x<-rank(x)
    if ((sum(w==1)==n)&(standardize)){
        vecX <- as.vector(x)
        extra <- max(xlabel) + 1
        na =  -93074815   # Consistency with mt.teststat
        nonpara<-"y"
        if (var.equal) test <- "t.equalvar"
        else
           test = "t"
        if (robust) test <- "wilcoxon"
        if (is.null(nrow(x))) traits <- 1
        else
            traits <- nrow(x)
        options <- c(test, "abs", "y")
        TestStat <- .C("get_stat_num_denum", as.double(vecX), as.integer(traits),
                 as.integer(n), as.integer(xlabel), as.double(na),
                 t.num = double(traits), t.denum = double(traits), as.character(options),
                 as.integer(extra), PACKAGE = "multtest")
         if (robust) {
             obs <- length(xlabel)
             lab1<-sort.int(unique.default(xlabel))[1]
             m=sum(xlabel==lab1)
             num <- TestStat$t.num*obs/(m*(obs-m))   # Conversion to same numerator value as original non C code function
             denom <- TestStat$t.denum*obs/(m*(obs-m)) # Conversion to same denominator value as original non C code function
         }
         else {
         num <- TestStat$t.num - psi0
         denom <- TestStat$t.denum
         }
     }
     else {
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
        if(sum(w)==1) df <- sum(w)
        else {
          df<- sum(w)-1
          }
        mm<-sum(w*x)/sum(w)
        if(standardize) denom <- sqrt((1/m+1/(n-m))*sum(w*(x-mm)^2)/df)
        else
          denom <- 1
			}
			else{
				if(standardize){
					df1<-sum(w1)-1
					df2<-sum(w2)-1
					df<-sum(w)-2
					if(var.equal) denom <- sqrt((1/m+1/(n-m))*(sum(w1*(sub1-m1)^2)+sum(w2*(sub2-m2)^2))/df)
					else
            denom <- sqrt((1/m*sum(w1*(sub1-m1)^2)/df1)+(1/(n-m)*sum(w2*(sub2-m2)^2)/df2))
				}
				else
					denom<-1
				}
			}
		}
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
		c(num,denom,snum)
	}
}

FX<-function(label,na.rm=TRUE,robust=FALSE){
	if(is.null(label))
		stop("A label variable is needed for this test")
	Samp<-1:length(label)
	function(x,w=NULL, samp=Samp){
		dep<-label[samp]
		if(is.null(w)){
			if(na.rm){
				drop<-is.na(x)|is.na(dep)
				x<-x[!drop]
				xlabel<-as.vector(dep[!drop])
			}
			else
				xlabel<-as.vector(dep)
    # Convert to 0,1,..k coding
    labs <- sort.int(unique.default(xlabel))
    num.levels <- length(labs)

    for (i in 1:num.levels){
        Index <- which(xlabel==labs[i])
        xlabel <- replace(xlabel, Index, i-1)
    }
    xlabel <- as.numeric(xlabel)
    # Check for at least 2 unique values in each group
    for (i in 1:length(labs)){
       if(length(unique.default(x[xlabel==labs[i]]))==1) stop("Only one unique value in bootstrap sample for one of the groups. Within group sum of squares is 0. This problem may be resolved if you try again with a different seed.")
    }
      if(robust) x<-rank(x)
    	n<-length(x)
      vecX <- as.vector(x)
      extra <- max(xlabel) + 1
      na =  -93074815 # Consistency with mt.teststat
      if (is.null(nrow(x))) traits <- 1
      else
        traits <- nrow(x)
      options <- c("f", "abs", "y")
      TestStat <- .C("get_stat_num_denum", as.double(vecX), as.integer(traits),
                 as.integer(n), as.integer(xlabel), as.double(na),
                 t.num = double(traits), t.denum = double(traits), as.character(options),
                 as.integer(extra), PACKAGE = "multtest")
      num <- TestStat$t.num
      denom <- TestStat$t.denum
    }
		else{
			if(length(w)!=length(x))
				stop("x and w must have same length")
			x[!is.finite(w)]<-NA
			if(na.rm){
				drop<-is.na(x)|is.na(dep)|is.na(w)
				x<-x[!drop]
				xlabel<-as.vector(dep[!drop])
				w<-w[!drop]
			}
			else
				xlabel<-as.vector(dep)
		# Convert to 0,1,..k coding
          labs <- sort.int(unique.default(xlabel))
          num.levels <- length(labs)

          for (i in 1:num.levels){
              Index <- which(xlabel==labs[i])
              xlabel <- replace(xlabel, Index, i-1)
          }
          xlabel <- as.numeric(xlabel)
          # Check for at least 2 unique values in each group
          for (i in 1:length(labs)){
             if(length(unique.default(x[xlabel==labs[i]]))==1) stop("Only one unique value in bootstrap sample for one of the groups. Within group sum of squares is 0. This problem may be resolved if you try again with a different seed.")
          }
			if(robust) x<-rank(x)
			n<-length(x)
      if(robust) x<-rank(x)
#TODO: how to deal with weights in F?
      vecX <- as.vector(x)
      extra <- max(xlabel) + 1
      na =  -93074815  # Consistency with mt.teststat
      if (is.null(nrow(x))) traits <- 1
      else
        traits <- nrow(x)
      options <- c("f", "abs", "y")
      TestStat <- .C("get_stat_num_denum", as.double(vecX), as.integer(traits),
                 as.integer(n), as.integer(xlabel), as.double(na),
                 t.num = double(traits), t.denum = double(traits), as.character(options),
                 as.integer(extra), PACKAGE = "multtest")
      num <- TestStat$t.num
      denom <- TestStat$t.denum
		}
		c(2*num,2*denom,1)
	}
}

#F statistic for block design with k treatments and l blocks
# One observation per block
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
# TODO: how to deal with weights in block f?
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
#
#F statistic for block design with k treatments and l blocks
# The observations are ordered by block, and within
#  each block, they are labeled using the integers 1 to k.
#Friedman statistic if robust=TRUE

twowayFX <-function(label,na.rm=TRUE,robust=FALSE){
	if(is.null(label))
		stop("A label variable is needed for this test")
	Samp<-1:length(label)
	function(x,w=NULL, samp=Samp){
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
      l <- length(gregexpr('12', paste(xlabel, collapse=""))[[1]])
			ublock<-1:l
      Breaks <- c(0,gregexpr(paste(c(k,1),collapse=""), paste(xlabel, collapse=""))[[1]], n)
      BlockNum <- sapply(1:l, function(x) Breaks[x+1]-Breaks[x])
      block <- unlist(sapply(1:l, function(x) rep(x,BlockNum[x])))
			if(robust){
			 	for(j in 1:l)
					x[block==j]<-rank(x[block==j])
			}
			m<-mean(x)
			mlab<-mblock<-mcell<-denom<-NULL
			for(i in 1:k){
				mlab[i]<-mean(x[xlabel==ulab[i]])
          for(j in 1:l){
					if(i==1)
						mblock[j]<-mean(x[block==ublock[j]])
            denom[(i-1)*l+j]<-sum((mean(x[xlabel==ulab[i] & block==ublock[j]])-mlab[i]-mblock[j]+m)^2)
				}
			}
      num<-sum(l*(mlab-m)^2)/(k-1)
			denom<-sum(denom)/((l-1)*(k-1))
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
      l <- length(gregexpr('12', paste(xlabel, collapse=""))[[1]])
			ublock<-1:l
      Breaks <- c(0,gregexpr(paste(c(k,1),collapse=""), paste(xlabel, collapse=""))[[1]], n)
      BlockNum <- sapply(1:l, function(x) Breaks[x+1]-Breaks[x])
      block <- unlist(sapply(1:l, function(x) rep(x,BlockNum[x])))
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
            denom[(i-1)*l+j]<-sum((mean(x[xlabel==ulab[i] & block==ublock[j]])-mlab[i]-mblock[j]+m)^2)
				}
			}
			num<-sum(l*(mlab-m)^2)/(k-1)
			denom<-sum(denom)/((l-1)*(k-1))
		}
		c(num,denom,1)
	}
}


## Z is a design *matrix*
## with variable of interest in first column
## and variables to adjust for in remaining columns
## Z is fixed for all columns of X
## gene expression is the outcome
lmX<-function(Z=NULL,n,psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	if(is.null(Z))
		Z<-matrix(1,n,1)
	else
		Z<-cbind(Z,rep(1,n))
	Samp<-1:n
	function(x,w=NULL, samp=Samp){
		covar<-Z[samp,]
		if(is.null(w)){
			if(na.rm){
				drop<-is.na(x)|rowSums(is.na(covar))
				covar<-covar[!drop,]
				x<-x[!drop]
			}
			covar<-as.matrix(covar)
			if(robust){
                          	autoload("rlm","MASS")
				out<-rlm(covar,x)
				out$df.residual<-length(x)-out$rank
			}
			else
				out<-lm.fit(covar,x)
			if(standardize) denom <- sqrt(sum(out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1])
         else denom <- 1
      if(denom==0) stop("Denominator of test statistic is 0 for a bootstrap sample. This problem may resuly from too small and sample size but may be resolved if you try again with a different seed.")
		}
		else{
			if(length(w)!=length(x))
				stop("x and w must have same length")
			x[!is.finite(w)]<-NA
			if(na.rm){
				drop<-is.na(x)|rowSums(is.na(covar))|is.na(w)
				covar<-covar[!drop,]
				x<-x[!drop]
				w<-w[!drop]
			}
			covar<-as.matrix(covar)
			if(robust){
                              	autoload("rlm","MASS")
				out<-rlm(covar,x,w)
				out$df.residual<-length(x)-out$rank
			}
			else
				out<-lm.wfit(covar,x,w)
			if(standardize) denom <- sqrt(sum(w*out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1])
         else denom <- 1
      if(denom==0) stop("Denominator of test statistic is 0 for a bootstrap sample. This problem may resuly from too small and sample size but may be resolved if you try again with a different seed.")
		}
		num<-out$coef[1]-psi0
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
		c(num,denom,snum)
	}
}


## gene expression is the covariate
## y is an outcome of interest
## Z is any intercept or other covariates
## Z changes for each row of X
lmY<-function(Y,Z=NULL,n,psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE){
	if(is.null(Y))
		stop("An outcome variable is needed for this test")
	if(length(Y)!=n)
		stop(paste("Dimension of outcome Y=",length(Y),", not equal dimension of data=",n,sep=""))
	if(is.null(Z))
		Z<-matrix(1,n,1)
	else
		Z<-cbind(Z,rep(1,n))
	Samp<-1:n
	function(x,w=NULL, samp=Samp){
		dep<-Y[samp]
		covar<-Z[samp,]
		covar<-cbind(x,covar)
		covar[!is.finite(w),]<-NA
		if(is.null(w)){
			if(na.rm){
				drop<-is.na(dep)|rowSums(is.na(covar))
				covar<-covar[!drop,]
				xy<-dep[!drop]
			}
			else
				xy<-dep
			if(robust){
                          	autoload("rlm","MASS")
                                out<-rlm(covar,xy)
				out$df.residual<-length(xy)-out$rank
			}
			else
				out<-lm.fit(covar,xy)

      if(standardize) denom <- sqrt(sum(out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1])
         else
           denom <- 1
      if(denom==0) stop("Denominator of test statistic is 0 for a bootstrap sample. This problem may resuly from too small and sample size but may be resolved if you try again with a different seed.")
		}
		else{
			if(length(w)!=length(x))
				stop("x and w must have same length")
			if(na.rm){
				drop<-is.na(dep)|rowSums(is.na(covar))|is.na(w)
				covar<-covar[!drop,]
				xy<-dep[!drop]
				w<-w[!drop]
			}
			else
				xy<-dep
			if(robust){
                            	autoload("rlm","MASS")
				out<-rlm(covar,xy,w)
				out$df.residual<-length(xy)-out$rank
			}
			else
				out<-lm.wfit(covar,xy,w)
			if(standardize) denom <- sqrt(sum(w*out$residuals^2)/out$df.residual)*sqrt(diag(chol2inv(out$qr$qr,size=out$rank))[1])
         else
            denom <- 1
      if(denom==0) stop("Denominator of test statistic is 0 for a bootstrap sample. This problem may resuly from too small and sample size but may be resolved if you try again with a different seed.")
		}
		num<-out$coef[1]-psi0
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
		c(num,denom,snum)
	}
}

## returns NA's if coxph fails
## strata is covariates to adjust for
coxY<-function(surv.obj,strata=NULL,psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",init=NULL,method="efron"){
	autoload("coxph","survival")
    	if(!inherits(surv.obj,"Surv"))  #covers NULL case
        	stop("Response must be a survival object")
	if(!is.null(strata))
		strat<-as.matrix(strata)
	else
		strat<-rep(1,nrow(surv.obj))
	strat<-strata(strat)
    	Samp<-1:nrow(surv.obj)
    	function(x,w=NULL, samp=Samp){
        	if(!is.null(w)&length(w)!=length(x))
            		stop("x and w must have same length")
		dep<-surv.obj[samp,]
		covar<-strat[samp]
		if(na.rm){
			drop<-is.na(x)
			if(!is.null(w))
				drop<-drop+is.na(w)
			if(!is.null(strat))
				drop<-drop+is.na(covar)
			x<-x[!drop]
			w<-w[!drop]
			covar<-covar[!drop]
			dep<-dep[!drop,]
		}
		if(sum(is.na(covar))){
			drop<-is.na(covar)
			x<-x[!drop]
			w<-w[!drop]
			covar<-covar[!drop]
			dep<-dep[!drop,]
		}
		design<-cbind(x,rep(1,length(x)))
		design[!is.finite(w),]<-NA
		control<-coxph.control()
		srvd<-try(coxph.fit(design,dep,strata=covar,init=init,control=control,weights=w,method=method,rownames=rownames(design)))
        	if(inherits(srvd,"try-error"))
	            return(c(NA,NA,NA))
	        if(standardize) denom <- sqrt(srvd$var[1,1])
             else denom <- 1
          if(denom==0) stop("Denominator of test statistic is 0 for a bootstrap sample. This problem may resuly from too small and sample size but may be resolved if you try again with a different seed.")
       	 	num<-srvd$coef[1]-psi0

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
		c(num,denom,snum)
	}
}

#function that applies stat.closure to (X,W)
get.Tn<-function(X,stat.closure,W=NULL){
	wapply(X,1,stat.closure,W)
}



