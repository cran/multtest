setClass("EBMTP",representation(statistic="numeric",
                              estimate="numeric",
                              sampsize="numeric",
                              rawp="numeric",
                              adjp="numeric",
                              reject="matrix",
                              rawdist="matrix",
                              nulldist="matrix",
                              nulldist.type="character",
                              marg.null="character",
                              marg.par="matrix",
                              label="numeric",
                              falsepos="matrix",
                              truepos="matrix",
                              errormat="matrix",  
                              EB.h0M="numeric",
                              prior="numeric",
                              prior.type="character",
                              lqv="numeric",
                              Hsets="matrix",
                              index="matrix",
                              call="call",
                              seed="integer"),
         prototype=list(statistic=vector("numeric",0),
         estimate=vector("numeric",0),
         sampsize=vector("numeric",0),
         rawp=vector("numeric",0),
         adjp=vector("numeric",0),
         reject=matrix(nr=0,nc=0),
         rawdist=matrix(nr=0,nc=0),
         nulldist=matrix(nr=0,nc=0),
         nulldist.type=vector("character",0),
         marg.null=vector("character",0),
         marg.par=matrix(nr=0,nc=0),
         label=vector("numeric",0),
         falsepos=matrix(nr=0,nc=0),
         truepos=matrix(nr=0,nc=0),
         errormat=matrix(nr=0,nc=0),
         EB.h0M=vector("numeric",0),
         prior=vector("numeric",0),
         prior.type=vector("character",0),
         lqv=vector("numeric",0),
         Hsets=matrix(nr=0,nc=0),
         index=matrix(nr=0,nc=0),
         call=NULL,
         seed=vector("integer",0)))

print.EBMTP<-function(x,...){
  call.list<-as.list(x@call)
  cat("\n")
  writeLines(strwrap("Multiple Testing Procedure",prefix="\t"))
  cat("\n")
  cat(paste("Object of class: ",class(x)))
  cat("\n")
  cat(paste("sample size =",x@sampsize,"\n"))
  cat(paste("number of hypotheses =",length(x@statistic),"\n"))
  cat("\n")
  cat(paste("test statistics =",ifelse(is.null(call.list$test),"t.twosamp.unequalvar",call.list$test),"\n"))
  cat(paste("type I error rate =",ifelse(is.null(call.list$typeone),"fwer",call.list$typeone),"\n"))
  nominal<-eval(call.list$alpha)
  if(is.null(eval(call.list$alpha))) nominal<-0.05
  cat("nominal level alpha = ")
  cat(nominal,"\n")
  cat(paste("multiple testing procedure =",ifelse(is.null(call.list$method),"common.cutoff",call.list$method),"\n"))
  cat("\n")
  cat("Call: ")
  print(x@call)
  cat("\n")
  cat("Slots: \n")
  snames<-slotNames(x)
  n<-length(snames)
  out<-matrix(nr=n,nc=4)
  dimnames(out)<-list(snames,c("Class","Mode","Length","Dimension"))
  for(s in snames) out[s,]<-c(class(slot(x,s)),mode(slot(x,s)),length(slot(x,s)),paste(dim(slot(x,s)),collapse=","))
  out<-data.frame(out)
  print(out)
  invisible(x)
}

### Put EBupdate last, since it is such a pain.
### Start with the rest of the other methods, and see what
### we want to keep/change from the MTP methods
# plot, EBMTP currently does not return cutoffs or confidence regions, so, leave 5 and 6 from MTP
# blank
if( !isGeneric("plot") ) setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod("plot","EBMTP",
	function(x,y="missing",which=1:4,caption=c("Rejections vs. Error Rate",
                                           "Ordered Adjusted p-values","Adjusted p-values vs. Statistics",
                                           "Unordered Adjusted p-values","Estimates & Confidence Regions",
                                           "Test Statistics & Cut-offs"),sub.caption = deparse(x@call,width.cutoff=500),
                   ask = prod(par("mfcol"))<length(which)&&dev.interactive(),
                   logscale=FALSE,top=10,...){
          call.list<-as.list(x@call)
          if(!inherits(x,"EBMTP")) stop("Use only with 'EBMTP' objects")
          if(is.null(which)) which<-1:4
          if(length(caption)==1) caption<-rep(caption,4)
          if(length(x@adjp)==0 & any(which)) stop("plot methods require adjusted p-values")
          #if(length(x@conf.reg)==0 & any(which==5)) stop("plot method 5 requires confidence regions")
          #if(length(x@cutoff)==0 & any(which==6)) stop("plot method 6 requires cut-offs")
          #go back to MTP method if we eventually want to put these in once cut-offs and conf reg
          #added into functionality, more for these below was deleted out. 
          if(!is.numeric(which) || any(which<1) || any(which>4)) stop("which must be in 1:4")
          show<-rep(FALSE,4)
          show[which]<-TRUE
          m<-length(x@adjp)
          if(top>m){
            warning("number of top hypotheses to plot exceeds total number of hypotheses - plotting less than requested number")
            top<-m
          }
	  ord<-order(x@adjp)
          if(any(show[2:4]) & logscale){
            pv<-(-log(x@adjp,10))
            pvlab<-"-log (base 10) Adjusted p-values"
          }
          else{
            pv<-x@adjp
            pvlab<-"Adjusted p-values"
          }
          one.fig<-prod(par("mfcol"))==1
          if(ask){
            op<-par(ask=TRUE)
            on.exit(par(op))
          }
          if(show[1]){
            nominal<-seq(0,1,by=0.05)
            r<-mt.reject(x@adjp,nominal)$r
            matplot(nominal,r,xlab="Type I error rate",
                    ylab="Number of rejected hypotheses",
                    type="l",...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[1],3,0.25)
          }
          if(show[2]){
            spval<-sort(pv)
            matplot(1:m,spval,xlab="Number of rejected hypotheses",
                    ylab=paste("Sorted",pvlab,sep=" "),type="l",...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[2],3,0.25)
          }
          if(show[3]){
            symb<-ifelse(length(pv)<100,"o",".")
            matplot(x@statistic,pv,xlab="Test statistics",
                    ylab=pvlab,type="p",pch=symb,...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[3],3,0.25)
          }
          if(show[4]){
            matplot(1:m,pv,xlab="Index",ylab=pvlab,type = "l", ...)
            if(one.fig) title(sub=sub.caption,cex.sub=0.5,...)
            mtext(caption[4],3,0.25)
          }
          if(!one.fig && par("oma")[3]>=1) mtext(sub.caption,outer=TRUE,cex=0.8)
          invisible()
          })


#summary
if( !isGeneric("summary") )
    setGeneric("summary", function(object, ...) standardGeneric("summary"))

setMethod("summary","EBMTP",
          function(object,...){
            call.list<-as.list(object@call)
            #cat(paste("EBMTP: ",ifelse(is.null(call.list$method),"common.cutoff",call.list$method),"\n"))
            cat("EBMTP: common.cutoff","\n") # always common.cutoff, even when being updated from MTP object
            err<-ifelse(is.null(call.list$typeone),"fwer",call.list$typeone)
            if(err=="gfwer") err<-paste(err," (k=",ifelse(is.null(call.list$k),0,call.list$k),")",sep="")
            if(err=="tppfp") err<-paste(err," (q=",ifelse(is.null(call.list$q),0.1,call.list$q),")",sep="")
	    cat(paste("Type I error rate: ",err,"\n"))
            cat(paste("prior: ",ifelse(is.null(call.list$prior),"conservative",call.list$prior),"\n\n"))
            nominal<-eval(call.list$alpha)
            if(is.null(nominal)) nominal<-0.05
            if(is.null(call.list$test)) test <- "t.twosamp.unequalvar"
            else test <- call.list$test
            if(test!="t.cor" & test!="z.cor") out1<-data.frame(Level=nominal,Rejections=apply(object@reject,2,sum),row.names=NULL)
            else{
              tmp <- rep(0,length(nominal))
              for(i in 1:length(nominal)) tmp[i] <- sum(object@adjp < nominal[i])
              out1 <- data.frame(Level=nominal,Rejections=tmp,row.names=NULL)
            }
            print(out1)
            cat("\n")
            out2<-get.index(object@adjp,object@rawp,abs(object@statistic))
            out3<-rn<-NULL
            if(!is.null(object@adjp)){
              out3<-rbind(out3,c(summary(object@adjp[!is.na(object@adjp)]),sum(is.na(object@adjp))))
              rn<-c(rn,"adjp")
            }
            if(!is.null(object@rawp)){
              out3<-rbind(out3,c(summary(object@rawp[!is.na(object@rawp)]),sum(is.na(object@rawp))))
              rn<-c(rn,"rawp")
            }
            if(!is.null(object@statistic)){
              out3<-rbind(out3,c(summary(object@statistic[!is.na(object@statistic)]),sum(is.na(object@statistic))))
            rn<-c(rn,"statistic")
            }
            if(!is.null(object@estimate)){
              out3<-rbind(out3,c(summary(object@estimate[!is.na(object@estimate)]),sum(is.na(object@estimate))))
              rn<-c(rn,"estimate")
            }
            rownames(out3)<-rn
            colnames(out3)[ncol(out3)]<-"NA's"
            print(out3)
            invisible(list(rejections=out1,index=out2,summaries=out3))
          })



if( !isGeneric("ebmtp2mtp") )
    setGeneric("ebmtp2mtp", function(object, ...) standardGeneric("ebmtp2mtp"))

setMethod("ebmtp2mtp","EBMTP",
          function(object,...){
            y<-new("MTP")
            slot(y,"statistic") <- object@statistic
            slot(y,"estimate") <- object@estimate
            slot(y,"sampsize") <- object@sampsize
            slot(y,"rawp") <- object@rawp
            slot(y,"adjp") <- object@adjp
            slot(y,"reject") <- object@reject
            slot(y,"rawdist") <- object@rawdist
            slot(y,"nulldist") <- object@nulldist
            slot(y,"nulldist.type") <- object@nulldist.type
            slot(y,"marg.null") <- object@marg.null
            slot(y,"marg.par") <- object@marg.par
            slot(y,"label") <- object@label
            slot(y,"index") <- object@index
            slot(y,"call") <- object@call
            slot(y,"seed") <- object@seed
            invisible(y)
          }
          )
            
setMethod("[","EBMTP",
          function(x,i,j=NULL,...,drop=FALSE){
            if(missing(i))
            i<-TRUE
            newx<-x
            slot(newx,"statistic")<-x@statistic[i]
            slot(newx,"estimate")<-x@estimate[i]
            slot(newx,"rawp")<-x@rawp[i]
            if(sum(length(x@adjp))) slot(newx,"adjp")<-x@adjp[i]
            if(sum(length(x@label))) slot(newx,"label")<-x@label[i]
            d<-dim(x@reject)
            dn<-dimnames(x@reject)
            if(sum(d)) slot(newx,"reject")<-matrix(x@reject[i,],nr=ifelse(i[1]==TRUE & !is.numeric(i),d[1],length(i)),nc=d[-1],dimnames=list(dn[[1]][i],dn[[2]]))
            if(sum(dim(x@nulldist))) slot(newx,"nulldist")<-x@nulldist[i,]
            if(sum(dim(x@marg.par))) slot(newx,"marg.par")<-x@marg.par[i,]
            if(sum(dim(x@rawdist))) slot(newx,"rawdist")<-x@rawdist[i,]
            if(sum(dim(x@falsepos))) slot(newx,"falsepos")<-x@falsepos[i,]
            if(sum(dim(x@truepos))) slot(newx,"truepos")<-x@truepos[i,]
            if(sum(dim(x@errormat))) slot(newx,"errormat")<-x@errormat[i,]
            slot(newx,"lqv")<-x@lqv[i]
            if(sum(dim(x@index))) slot(newx,"index")<-x@index[i,]
	    invisible(newx)
          })

setMethod("as.list","EBMTP",
          function(x,...){
            snames<-slotNames(x)
            n<-length(snames)
            lobj<-list()
            for(i in 1:n) lobj[[i]]<-slot(x,snames[i])
            names(lobj)<-snames
            invisible(lobj)
          })


if( !isGeneric("EBupdate") )
    setGeneric("EBupdate", function(object, ...) standardGeneric("EBupdate"))

setMethod("EBupdate","EBMTP",
          function(object,formula.="missing",alternative="two.sided",typeone="fwer",
          k=0,q=0.1,alpha=0.05,smooth.null=FALSE,
          method="common.cutoff",prior="conservative",bw="nrd",kernel="gaussian",
          get.adjp=TRUE,nulldist="boot.cs",keep.rawdist=FALSE,keep.nulldist=TRUE,
          keep.falsepos=FALSE,keep.truepos=FALSE,keep.errormat=FALSE,keep.Hsets=FALSE,
          marg.null=object@marg.null,marg.par=object@marg.par,ncp=NULL,
          keep.label=TRUE,...,evaluate=TRUE){
            p <- length(object@statistic)
            m <- length(object@statistic)
            B <- dim(object@nulldist)[2]
            if(sum(object@rawdist)!=0) B <- dim(object@rawdist)[2]
            ## checking
            #Error rate
            ERROR<-c("fwer","gfwer","tppfp","fdr")
            typeone<-ERROR[pmatch(typeone,ERROR)]
            if(is.na(typeone)) stop(paste("Invalid typeone, try one of ",ERROR,sep=""))
            if(any(alpha<0) | any(alpha>1)) stop("Nominal level alpha must be between 0 and 1.")
            nalpha<-length(alpha)
            reject<-
              if(nalpha) array(dim=c(p,nalpha),dimnames=list(names(object@rawp),paste("alpha=",alpha,sep="")))
              else matrix(nr=0,nc=0)

            if(typeone=="fwer"){
              if(length(k)>1) k<-k[1]
              if(sum(k)!=0) stop("FWER control, by definition, requires k=0.  To control k false positives, please select typeone='gfwer'.")
            }
     
            if(typeone=="gfwer"){
              if(length(k)>1){
                k<-k[1]
                warning("Can only compute gfwer adjp for one value of k at a time (using first value). Use EBupdate() to get results for other values of k.")
              }
              if(k<0) stop("Number of false positives can not be negative.")
              if(k>=p) stop(paste("Number of false positives must be less than number of tests=",p,sep=""))
            }

            if(typeone=="tppfp"){
              if(length(q)>1){
                q<-q[1]
                warning("Can only compute tppfp adjp for one value of q at a time (using first value). Use EBupdate() to get results for other values of q.")
              }
              if(q<0) stop("Proportion of false positives, q, can not be negative.")
              if(q>1) stop("Proportion of false positives, q, must be less than 1.")
            }
            
            #methods
            METHODS<-c("common.cutoff","common.quantile")
            method<-METHODS[pmatch(method,METHODS)]
            if(is.na(method)) stop(paste("Invalid method, try one of ",METHODS,sep=""))
            if(method=="common.quantile") stop("Common quantile procedure not currently implemented.  Common cutoff is pretty good, though.")
            
            #prior
            PRIORS<-c("conservative","ABH","EBLQV")
            prior<-PRIORS[pmatch(prior,PRIORS)]
            if(is.na(prior)) stop(paste("Invalid prior, try one of ",PRIORS,sep=""))

            #get args from previous call
            call.list<-as.list(object@call)

            if(is.null(call.list$test)) test<-"t.twosamp.unequalvar" #default
            else test<-call.list$test
            ### nulldistn
            ### Preserve the old null dist, if kept (i.e., could have alternatively kept raw dist)
            nulldistn <- object@nulldist
            if(object@nulldist.type=="perm") stop("No way to update objects which originally used the permutation distribution. No available options for storing nulldist.  Rawdist can only be stored for bootstrap distribution.")
            ### For boot.qt, make sure values of marg.null and marg.par, if set previously, are kept.
            ### Otherwise, these become null, but the original values are set here before proceeding.
            prev.marg.null <- object@marg.null
            prev.marg.par <- object@marg.par

            if(!ncol(object@nulldist) & !ncol(object@rawdist)) stop("Update method requires either keep.raw and/or keep.null=TRUE in original call to MTP")
            nulldist<- # just setting character value of what nulldist should be
               if(is.null(call.list$nulldist)) "boot.cs"
               else call.list$nulldist

         ## new call
               newcall.list<-as.list(match.call())
               changed<-names(call.list)[names(call.list)%in%names(newcall.list)]
               changed<-changed[changed!=""]
               added<-names(newcall.list)[!(names(newcall.list)%in%names(call.list))]
               added<-added[added!="x"]
               for(n in changed) call.list[[n]]<-newcall.list[[n]]
               for(n in added) call.list[[n]]<-newcall.list[[n]]
               newcall<-as.call(call.list)
               ### NB can still use "call.list" to help with what has been changed.
               df <- marg.par
               call.list$marg.par <- df
               
         ## return call if evaluate is false
               if(!evaluate) return(newcall)

         ## else redo MTP
            else{
              num<-object@estimate
              snum<-1
              if(alternative=="two.sided"){
                snum<-sign(num)
                num<-abs(num)
              }
              if(alternative=="less"){
                snum<-(-1)
                num<-(-num)
              }

              if(object@nulldist.type!="boot.qt"){
                marg.null = vector("character",length=0)
                marg.par = matrix(nr=0,nc=0)
              }
              if("alternative" %in% changed | "alternative" %in% added) alternative <- call.list$alternative
              if("marg.null" %in% changed | "marg.null" %in% added) marg.null <- call.list$marg.null
              if("marg.par" %in% changed | "marg.par" %in% added){
                  marg.par <- call.list$marg.par
                  if(is.numeric(marg.par) & !is.matrix(marg.par)) marg.par <- matrix(rep(marg.par,length(object@statistic)),nr=length(object@statistic),nc=length(marg.par),byrow=TRUE)
                }
              if("perm.mat" %in% changed | "perm.mat" %in% added) perm.mat <- call.list$perm.mat
              if("ncp" %in% changed | "ncp" %in% added) ncp <- call.list$ncp
              if("MVN.method" %in% changed | "MVN.method" %in% added | "penalty" %in% changed | "penalty" %in% added |"ic.quant.trans" %in% changed | "ic.quant.trans" %in% added) stop("Changing 'MVN.method', 'ic.quant.trans' or 'penalty' requires new calculation of null distribution using nulldist='ic'.  Please use a new call to EBMTP.")
         ### Check value of nulldist in this case
              if("nulldist" %in% changed | "nulldist" %in% added) {
                nulldist <- call.list$nulldist
         ### Otherwise, nulldist keeps the old/default value in the original call.list, not the updated one.
                if(nulldist=="perm") stop("Calls to update() cannot include changes involving the permutation distribution. Please try a separate call to MTP() with nulldist='perm'")
                if(object@nulldist.type=="ic") stop("You cannot update an influence curve null distribution to another choice of null distribution.  Valid only for changes in the bootstrap distribution when keep.rawdist=TRUE.  Please try a separate call to MTP() if nulldist='boot' or 'perm' desired. Changing 'MVN.method', 'ic.quant.trans' or 'penalty' also requires new calculation of null distribution using nulldist='ic'")
                if(nulldist=="ic") stop("Calls to update() cannot include changes involving the influence curve null distribution. Please try a separate call to MTP() with nulldist='ic'")
                if(!ncol(object@rawdist)) stop("Calls to update() involving changes in bootstrap-based null distributions require keep.rawdist=TRUE")

    ### Just recompute (bootstrap-based) nulldistn - way easier this way (with keep.raw=TRUE)
    ### "Easy" ones first.  Need to get tau0 and theta0.
              if(nulldist=="ic"){
                marg.null = vector("character",length=0)
                marg.par = matrix(nr=0,nc=0)
              }
                
              if(nulldist=="boot" | nulldist=="boot.cs" | nulldist=="boot.ctr"){
                marg.null = vector("character",length=0)
                marg.par = matrix(nr=0,nc=0)
                tau0<-1
                theta0<-0
                if(test=="f"){
                  theta0<-1
                  tau0<-2/(length(unique(object@label))-1)
                }
                if(test=="f.twoway"){
                  theta0<-1
                  tau0 <- 2/((length(unique(object@label))*length(gregexpr('12', paste(object@label, collapse=""))[[1]]))-1)
                }
                if(nulldist=="boot") nulldistn <- center.scale(object@rawdist, theta0, tau0, alternative)
                if(nulldist=="boot.cs") nulldistn <- center.scale(object@rawdist, theta0, tau0, alternative)
                if(nulldist=="boot.ctr") nulldistn <- center.only(object@rawdist, theta0, alternative)
              }

              if(nulldist=="boot.qt"){
                if("marg.null" %in% changed | "marg.null" %in% added) marg.null <- call.list$marg.null
                else marg.null <- NULL
                if("marg.par" %in% changed | "marg.par" %in% added){
                  marg.par <- call.list$marg.par
                  if(is.numeric(marg.par) & !is.matrix(marg.par)) marg.par <- matrix(rep(marg.par,length(object@statistic)),nr=length(object@statistic),nc=length(marg.par),byrow=TRUE)
                }
                else marg.par <- NULL
      
        ### If these additional args are changed or added, these will be the new defaults, but they will not be NULL
                ### Cannot be NULL for object defn.
                ncp <- if(is.null(call.list$ncp)) 0
                perm.mat <- if(is.null(call.list$perm.mat)) NULL
                if(!is.null(perm.mat)){
                  if(length(object@statistic)!=dim(perm.mat)[1]){ stop("Permutation and bootstrap matrices must have same number of rows (hypotheses).")
                                                                }
                }

                nstats <- c("t.twosamp.unequalvar","z.cor","lm.XvsZ","lm.YvsXZ","coxph.lmYvsXZ")
                tstats <- c("t.onesamp","t.twosamp.equalvar","t.pair","t.cor")
                fstats <- c("f","f.block","f.twoway")
         # If default (=NULL), set values of marg.null to pass on.
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
                          t.onesamp = object@sampsize-1,
                          t.twosamp.equalvar = object@sampsize-2,
                          t.twosamp.unequalvar = c(0,1),
                          t.pair = object@sampsize-2,
                          f = c(length(is.finite(unique(object@label)))-1,object@sampsize-length(is.finite(unique(object@label)))),
                          f.twoway = {
                            c(length(is.finite(unique(object@label)))-1,object@sampsize-(length(is.finite(unique(object@label)))*length(gregexpr('12', paste(y, collapse=""))[[1]]))-2)
                            },
                          lm.XvsZ = c(0,1),
                          lm.YvsXZ = c(0,1),
                          coxph.YvsXZ = c(0,1),
                          t.cor = object@sampsize-2,
                          z.cor = c(0,1)
                          )
                  marg.par <- matrix(rep(marg.par,length(object@statistic)),nr=length(object@statistic),nc=length(marg.par),byrow=TRUE)
        }
                else{ # Check that user-supplied values of marg.par make sense (marg.par != NULL)
                  if((marg.null=="t" | marg.null=="f") & any(marg.par[,1]==0)) stop("Cannot have zero df with t or F distributions. Check marg.par settings")
                  if(marg.null=="t" & dim(marg.par)[2]>1) stop("Too many parameters for t distribution.  marg.par should have length 1.")
                  if((marg.null=="f" | marg.null=="normal") & dim(marg.par)[2]!=2) stop("Incorrect number of parameters defining marginal null distribution.  marg.par should have length 2.")
                }
                nulldistn <- quant.trans(object@rawdist, marg.null, marg.par, ncp, alternative, perm.mat)
              }
              }

              ### Cool. Now pick up where we left off.
              ##performing multiple testing
              #rawp values
              obs<-rbind(num,object@estimate/object@statistic,sign(object@estimate))
              rawp<-apply((obs[1,]/obs[2,])<=nulldistn,1,mean)
		     if(smooth.null & min(rawp,na.rm=TRUE)==0){
                       zeros<-rawp==0
                       if(sum(zeros)==1){
                         den<-density(nulldistn[zeros,],to=max(obs[1,zeros]/obs[2,zeros],nulldistn[zeros,],na.rm=TRUE),na.rm=TRUE)
                         rawp[zeros]<-sum(den$y[den$x>=(obs[1,zeros]/obs[2,zeros])])/sum(den$y)
                       }
                       else{
                         den<-apply(nulldistn[zeros,],1,density,to=max(obs[1,zeros]/obs[2,zeros],nulldistn[zeros,],na.rm=TRUE),na.rm=TRUE)
                         newp<-NULL
                         stats<-obs[1,zeros]/obs[2,zeros]
                         for(i in 1:length(den)) newp[i]<-sum(den[[i]]$y[den[[i]]$x>=stats[i]])/sum(den[[i]]$y)
                         rawp[zeros]<-newp
                       }
                       rawp[rawp<0]<-0
                     }

              #c, cr, adjp - this is where the function gets a lot different from MTP.
              ### Begin nuts and bolts of EB here.
t
              ### Set G function of type I error rates
              error.closure <- switch(typeone, fwer=G.VS(V,S=NULL,tp=TRUE,bound=0),
                                      gfwer=G.VS(V,S=NULL,tp=TRUE,bound=k),
                                      tppfp=G.VS(V,S,tp=TRUE,bound=q),
                                      fdr=G.VS(V,S,tp=FALSE,bound=NULL)
                                      )

              ### Generate guessed sets of true null hypotheses
              ### This function relates null and full densities.  Sidedness should be accounted for above.
              statistic <- (obs[3,]*obs[1,]/obs[2,]) #observed, with sign
              Tn <- obs[1,]/obs[2,]  # for sidedness, matching with mulldistn
              
              H0.sets <- Hsets(Tn, nullmat=nulldistn, bw, kernel, prior=prior, B=dim(object@nulldist)[2], rawp=object@rawp) 
              EB.h0M <- H0.sets$EB.h0M
              prior.type <- prior
              prior.val <- H0.sets$prior
              lqv <- H0.sets$pn.out
              H0.sets <- H0.sets$Hsets.mat
              
              m <- length(Tn)
              ### B defined in global environment
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
              Vn <- object
              Vn <- .Call(VScount,as.numeric(Z.nulls),as.numeric(cutoffs),as.integer(m),
                as.integer(B),as.integer(clen),NAOK=TRUE)
              cat("\n")
              Vn <- matrix(Vn, nr=clen, nc=B)
              
              if(typeone=="fwer" | typeone=="gfwer") Sn <- NULL
              else{
                cat("counting guessed true positives...", "\n")
                Sn <- .Call(VScount,as.numeric(Tn.mat),as.numeric(cutoffs),as.integer(m),
                  as.integer(B),as.integer(clen),NOAK=TRUE)
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
              if(keep.errormat) G <- G[rev.order,]
              else G <- matrix(0,nr=0,nc=0)
              if(!keep.Hsets) H0.sets <- matrix(0,nr=0,nc=0)
              
              # No confidence regions, but vector of rejections logicals, and cutoff, if applicable
              ### Generate matrix of rejection logicals.
              EB.reject <- matrix(rep(0,m),nr=m,nc=length(alpha))
              dimnames(EB.reject) <- list(rownames(object@nulldist),paste("alpha", alpha, sep=""))
              if(nalpha) for(a in 1:nalpha) EB.reject[,a]<-adjp<=alpha[a]
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
              if(!keep.nulldist) nulldistn <-matrix(nr=0,nc=0)
              if(keep.rawdist==FALSE) object@rawdist<-matrix(nr=0,nc=0)
                out<-new("EBMTP",statistic=object@statistic,estimate=object@estimate,
                sampsize=object@sampsize,rawp=rawp,adjp=adjp,
                reject=EB.reject,rawdist=object@rawdist,nulldist=nulldistn,
                nulldist.type=nulldist,marg.null=marg.null,marg.par=marg.par,label=object@label,
                falsepos=Vn,truepos=Sn,errormat=G,Hsets=H0.sets,EB.h0M=EB.h0M,
                prior=prior.val,prior.type=prior.type,lqv=lqv,
                index=object@index,call=newcall,seed=object@seed)
		return(out)
               } #re else redo MTP
             } # re function
             ) # re set method
