.mt.BLIM<-2^30
.mt.naNUM<- -93074815
.mt.RandSeed<-3455660
#the maxim number of setting of the permutation, it's not resettable
#in the current version. the numer comes from the largest
#integer can be 2^32, while we need to exlcude one sign bit, and
#to exclude another bit for safety.


#dyn.load("multtest.so")
#X is a matrix data
#classlabel is a vector
mt.teststat<-function(X,classlabel,test="t",na=.mt.naNUM,nonpara="n")
{
    if(is.factor(classlabel)) classlabel<-unclass(classlabel)-1
    extra<-max(classlabel)+1
    mt.checkothers(na=na,nonpara=nonpara)
    tmp<-mt.transformX(X,classlabel,test,na,nonpara)
    options<-c(test,"abs","y"); #"abs"  and "y" has no meaning here
    res<-.C("get_stat",as.double(tmp$X),as.integer(tmp$m),
               as.integer(tmp$n),as.integer(tmp$classlabel),as.double(na),
               teststat=double(tmp$m),as.character(options),
               as.integer(extra), PACKAGE="multtest")$teststat
    res[abs(res)>=0.9*1e20]<-NA
    res
}
mt.teststat.num.denum<-function(X,classlabel,test="t",na=.mt.naNUM,nonpara="n")
{
    extra<-max(classlabel)+1
    mt.checkothers(na=na,nonpara=nonpara)
    tmp<-mt.transformX(X,classlabel,test,na,nonpara)
    options<-c(test,"abs","y"); #"abs"  and "y" has no meaning here
    teststat<-.C("get_stat_num_denum",as.double(tmp$X),as.integer(tmp$m),
	       as.integer(tmp$n),as.integer(tmp$classlabel),as.double(na),
	       t.num=double(tmp$m),t.denum=double(tmp$m),as.character(options),
               as.integer(extra), PACKAGE="multtest")

    res<-cbind(teststat.num=teststat$t.num,teststat.denum=teststat$t.denum)
    mt.niceres(res,X)
}
  mt.maxT<-function(X,classlabel,test="t",side="abs",
                  fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
{
    if(is.factor(classlabel)) classlabel<-unclass(classlabel)-1
    extra<-max(classlabel)+1
    mt.checkothers(side=side,fixed.seed.sampling=fixed.seed.sampling,B=B,na=na,nonpara=nonpara)
    tmp<-mt.transformX(X,classlabel,test,na,nonpara)
    newB<-mt.getmaxB(classlabel,test,B)
    if(B==0||newB<B)
      fixed.seed.sampling<-"n" #as we're doing complete premutation
    options<-c(test,side,fixed.seed.sampling);
    res<-.C("get_maxT",as.double(tmp$X),as.integer(tmp$m),
	    as.integer(tmp$n),as.integer(tmp$classlabel),as.double(na),
	    t=double(tmp$m),p=double(tmp$m),adjP=double(tmp$m),
	    as.integer(newB),index=integer(tmp$m),as.character(options),
               as.integer(extra), PACKAGE="multtest")

    res<-cbind(index=res$index,teststat=res$t,rawp=res$p,adjp=res$adjP)
    mt.niceres(res,X,res[,1])
}
mt.minP<-function(X,classlabel,test="t",side="abs",
                  fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
{
    if(is.factor(classlabel)) classlabel<-unclass(classlabel)-1
    extra<-max(classlabel)+1
    mt.checkothers(side=side,fixed.seed.sampling=fixed.seed.sampling,B=B,na=na,nonpara=nonpara)
    tmp<-mt.transformX(X,classlabel,test,na,nonpara)
    newB<-mt.getmaxB(classlabel,test,B)
    if(B==0||newB<B)
      fixed.seed.sampling<-"n" #as we're doing complete premutation
    options<-c(test,side,fixed.seed.sampling);
    res<-.C("get_minP",as.double(tmp$X),as.integer(tmp$m),
	    as.integer(tmp$n),as.integer(tmp$classlabel),as.double(na),
	    t=double(tmp$m),p=double(tmp$m),adjP=double(tmp$m),
            plower=double(tmp$m),as.integer(newB),index=integer(tmp$m),
            as.character(options),as.integer(extra), PACKAGE="multtest")

    res<-cbind(index=res$index,teststat=res$t,rawp=res$p,adjp=res$adjP,plower=res$plower)
    mt.niceres(res,X,res[,1])
}
mt.sample.teststat<-function(V,classlabel,test="t",fixed.seed.sampling="y",
                       B=10000,na=.mt.naNUM,nonpara="n")
{
  extra<-max(classlabel)+1
  mt.checkothers(fixed.seed.sampling=fixed.seed.sampling,B=B,na=na,nonpara=nonpara)
  tmp<-mt.transformV(V,classlabel,test,na,nonpara)
  newB<-mt.getmaxB(classlabel,test,B)
  if(B==0||newB<B)
    fixed.seed.sampling<-"n" #as we're doing complete premutation
  options<-c(test,"abs",fixed.seed.sampling);#the "abs" has no meaing here.
  res<-t(.C("get_samples_T",as.double(tmp$V),as.integer(tmp$n),
          as.integer(tmp$classlabel),T=double(newB),as.double(na),
          as.integer(newB),as.character(options),as.integer(extra),
            PACKAGE="multtest")$T)
  res[abs(res)>=0.9*1e20]<-NA
  res
}
mt.sample.rawp<-function(V,classlabel,test="t",side="abs",
                       fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
{
  extra<-max(classlabel)+1
  mt.checkothers(side=side,fixed.seed.sampling=fixed.seed.sampling,B=B,na=na,nonpara=nonpara)
  tmp<-mt.transformV(V,classlabel,test,na,nonpara)
  newB<-mt.getmaxB(classlabel,test,B)
  if(B==0||newB<B)
    fixed.seed.sampling<-"n" #as we're doing complete premutation
  options<-c(test,side,fixed.seed.sampling);
  res<-.C("get_samples_P",as.double(tmp$V),as.integer(tmp$n),
          as.integer(tmp$classlabel), P=double(newB),as.double(na),
          as.integer(newB),
          as.character(options),as.integer(extra), PACKAGE="multtest")$P
  res[abs(res)>=0.9*1e20]<-NA
  res
}
mt.sample.label<-function(classlabel,test="t",
                            fixed.seed.sampling="y",B=10000)
{
  extra<-max(classlabel)+1
  tmp<-mt.transformL(classlabel,test)
  mt.checkothers(fixed.seed.sampling=fixed.seed.sampling,B=B)
  newB<-mt.getmaxB(classlabel,test,B)
  if(B==0||newB<B)
    fixed.seed.sampling<-"n" #as we're doing complete premutation
  options<-c(test,"abs",fixed.seed.sampling); #the "abs" has no meaing here
  res<-.C("get_sample_labels",as.integer(tmp$n),as.integer(tmp$classlabel),
          as.integer(newB), S=integer(tmp$n*newB),as.character(options),
          as.integer(extra), PACKAGE="multtest")$S
  resl<-matrix(res,nrow=tmp$n)
  if(test=="pairt"){
    #restore the original classlabelling
    resn<-matrix(0,nrow=2*tmp$n,ncol=newB)
    for(i in c(1:tmp$n))
      for(j in c(1:newB)){
        if(resl[i,j])
          resn[2*i,j]<-1
        else resn[2*i-1,j]<-1
      }
    resl<-resn
  }
  t(resl)
}
#used for private function
mt.checkclasslabel<-function(classlabel,test)
{
  classlabel<-as.integer(classlabel)
  if((!is.character(test))||(!is.vector(test))||(length(test)>1)
      ||(!any(test==c("t","f","blockf","pairt","wilcoxon","t.equalvar"))))
     stop(paste("your setting of test is",test,"\nthe test needs to be a single character from c('t',f','blockf','pairt','wilcoxon','t.equalvar')"))
  if((!is.integer(as.integer(classlabel))) ||(!is.vector(classlabel)))
     stop("classlabel needs to be just a vector of integers")
  if(any(test==c("t","wilcoxon","t.equalvar"))){
    x<-sum(classlabel==0)
    y<-sum(classlabel==1)
    if((x==0)||(y==0)||(x+y<length(classlabel)))
      stop(paste("in t test, every number in class label needs to be 0 or 1 and neither of the 0 set or 1 set can be empty set\n",
                 "The folllowing is your setting of classlabel",classlabel,"\n"))
  }
  if(test=="f"){
      tab <- table(classlabel)
      tab <- tab[tab>0]
      if(length(tab)<2)
          stop(paste("in F test, we need at least two groups\n",
                     "Your setting of classlabel is", classlabel,
                   "\n"))
      if(sum(tab)-length(tab)<2)
          stop(paste("Insufficient df for denominator of F",
                     "the settings are", classlabel, "\n"))
  }
  if(test=="pairt"){
    K<-max(classlabel)
    if(K!=1)
      stop(paste("in paired t test, we only handle two groups\n",
                 "your classlabel=",classlabel,"\n"))
    if(length(classlabel)%%2==1)
      stop(paste("the classlabel length must be an even number in the paired t\n","your classlabel=",classlabel,"\n"))
    halfn<-length(classlabel)%/%2
    for(i in c(1:halfn)){
      cur<-classlabel[(2*i-1):(2*i)]
      if((sum(cur==0)==0)||(sum(cur==1)==0))
        stop(paste("Some errors in specifying classlabel for the paired t test for the block",i,"located at","(",2*i-1,2*i,")\n",
                   "your classlabel=",classlabel,"\n"))

    }
  }
  if(test=="blockf"){
    K<-max(classlabel)
    if(K<1)
      stop(paste("in blockF test, we need at least two groups\n",
                 "your classlabel=",classlabel,"\n"))
    if(length(classlabel)%%(K+1)>0)
      stop(paste("the classlabel length must be the multiple of the number of treatments in the block test\n","your classlabel=",classlabel,"\n"))
     B<-length(classlabel)%/%(K+1)
    for(i in c(1:B)){
      cur<-classlabel[c((K+1)*(i-1)+1):((K+1)*i)]
      #to check if cur is a permutation of c(0,1,..,K)
      for(j in c(0:K))
        if(sum(cur==j)==0)
          stop(paste("the classlabel has some errors for the blockf test at block",i,"located at",
               "(",(K+1)*(i-1)+1,(K+1)*i,")","There is no elements =",j,"within this block\n","your classlabel=",classlabel,"\n"))
    }
  }
}
mt.checkX<-function(X,classlabel,test){
  if((!is.matrix(X)) || !(is.numeric(X)))
     stop(paste("X needs to be a matrix\n","your X=",X,"\n"))
  if(ncol(X)!=length(classlabel))
    stop(paste("the number of column of X needs to be the same as the lengtho of classlabel\n","your X=",X,"\n your classlabel is",classlabel,"\n"))
  mt.checkclasslabel(classlabel,test)
}
mt.checkV<-function(V,classlabel,test){
  if((!is.vector(V)) || !(is.numeric(V)))
     stop(paste("V needs to be a vector\n","your V=",V,"\n"))
  if(length(V)!=length(classlabel))
    stop("the length of V needs to be the same as the length of classlabel\n",
         "your V=",V,"\n your classlabel=",classlabel,"\n")
  mt.checkclasslabel(classlabel,test)
}

mt.checkothers<-function(side="abs",fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
{
  if((length(B)>1) || !(is.integer(as.integer(B))) ||(!is.vector(B)))
     stop(paste("B needs to be just a integer\n","your B=",B,"\n"))
   if(B<0)
     stop(paste("the number of Permutations (B) needs to be positive\n, If you want to complete permutation, just specify B as any number greater than the maximum number of permutation\n","your B=",B))
  if((length(na)>1) || !(is.numeric(na)) ||(!is.vector(na)))
     stop(paste("na needs to be just a number\n","your na=",na,"\n"))
  if((!is.character(side))||(!is.vector(side))||(length(side)>1)
     ||(!any(side==c("upper","abs","lower"))))
    stop(paste("the side needs to be a single character from c('upper','abs','lower')\n","your side=",side,"\n"))
  if((!is.character(fixed.seed.sampling))||(!is.vector(fixed.seed.sampling))||(length(fixed.seed.sampling)>1)
     ||(!any(fixed.seed.sampling==c("y","n"))))
    stop(paste("the fixed.seed.sampling needs to be a single character from c('y','n')\n","your fixed.sampling=",fixed.seed.sampling,"\n"))
  if((!is.character(nonpara))||(!is.vector(nonpara))||(length(nonpara)>1)
     ||(!any(nonpara==c("y","n"))))
    stop(paste("the nonpara needs to be a single character from c('y','n')\n","your nonpara=",nonpara,"\n"))
}
mt.transformX<-function(X,classlabel,test,na,nonpara)
{
  X<-mt.number2na(data.matrix(X),na)
  mt.checkX(X,classlabel,test)
  n<-ncol(X)
  if(test=="pairt"){
    if(n%%2==1)
      stop(paste("the number of columns for X must be an even number in the paired t  test\n","your X=",X,"\n your classlabel=",classlabel,
                 "\n your test=",test,"\n"))
    halfn<-n%/%2;
    evendata<-X[,c(1:halfn)*2]
    odddata<-X[,c(1:halfn)*2-1]
    vecX<-(evendata-odddata)
    vecX<-data.matrix(vecX)
  }else{
    vecX<-data.matrix(X)
  }
  if(test=="wilcoxon"||nonpara=="y"){
    for(i in c(1:nrow(vecX))){
      vecX[i,]<-rank(vecX[i,])
    }
  }
  vecX<-mt.na2number(c(vecX),na)
  newL<-mt.transformL(classlabel,test)
  list(X=vecX,m=nrow(X),n=newL$n,classlabel=newL$classlabel)
}
mt.transformV<-function(V,classlabel,test,na,nonpara)
{
  V<-mt.number2na(as.double(V),na)
  mt.checkV(V,classlabel,test)
  n<-length(classlabel)
  if(test=="pairt"){
    if(n%%2==1)
      stop(paste("the number of columns for V must be an even number in the paired t test\n","your V=",V,"\n your classlabel=",classlabel,
                 "\n your test=",test,"\n"))
    halfn<-n%/%2
    evendata<-V[c(1:halfn)*2]
    odddata<-V[c(1:halfn)*2-1]
    newV<-c(evendata-odddata)
  }
  else{
    newV<-V
  }
  if(test=="wilcoxon"||nonpara=="y"){
    newV<-rank(newV)
    }
  newL<-mt.transformL(classlabel,test)
  list(V=mt.na2number(newV,na),n=newL$n,classlabel=newL$classlabel)
}
mt.transformL<-function(classlabel,test)
{
  classlabel<-as.integer(classlabel)
  mt.checkclasslabel(classlabel,test)
  n<-length(classlabel)
  newL<-classlabel
  if(test=="pairt"){
    if(n%%2==1)
      stop(paste("the length of classlabel must be an even number in the pair t\n","your classlabel=",classlabel,"\n your test=",test="\n"))
    halfn<-n%/%2;
    n<-halfn
    newL<-rep(0,n);
    for(i in c(1:n)){
      newL[i]<-classlabel[2*i]
    }
  }
  list(classlabel=newL,n=n)
}
#this functions finds the maximum number of permutation
#if the the initial B=0, or initial B greater than the maximum number of
#permutation maxB, it will return all possible of number of permutation.
mt.getmaxB<-function(classlabel,test,B, verbose=FALSE)
{
  if(B>.mt.BLIM)
    stop(paste("The setting of B=",B,"is too large, Please set B<",.mt.BLIM,"\n"))
  n<-length(classlabel)
  if(test=="pairt"){
    maxB<-2^(n%/%2)
  }
  if(any(test==c("t","f","wilcoxon","t.equalvar"))){
    k<-max(classlabel)
    maxB<-1
    curn<-n
    for(i in c(0:k)){
      nk<-sum(classlabel==i)
      for(j in c(1:nk)){
        maxB<-maxB*curn/j
        curn<-curn-1
      }
    }
  }
  if(test=="blockf"){
    k<-max(classlabel)
    maxB<-1
    for(i in c(1:(k+1))){
      maxB<-maxB*i
    }
    maxB<-maxB^(n%/%(k+1))
  }
  #finished the computing of maxB
  if((B==0)&(maxB>.mt.BLIM)){
    stop(paste("The complete enumeration is too big",maxB,
               "is too large, Please set random permutation\n"))
  }
  if((B>maxB)||(B==0)){
    if(verbose) cat("We'll do complete enumerations\n")
    return(maxB)
  }
  return(B)
}

mt.na2number<-function(x,na){
  y<-x
  y[is.na(y)]<-na
  y
}
mt.number2na<-function(x,na){
  y<-x
  y[y==na]<-NA
  y
}
#patched from the new version
mt.niceres<-function(res,X,index){
  newres<-res
  name<-rownames(X,do.NULL=FALSE,prefix="")
  if(missing(index)) {
    rownames(newres)<-name
  }else {
    rownames(newres)<-name[index]
  }
  newres[abs(newres)>=0.9*1e20]<-NA
  data.frame(newres)
}







