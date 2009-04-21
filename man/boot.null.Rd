\name{boot.null}
\alias{boot.null}
\alias{boot.resample}
\alias{center.scale}
\alias{center.only}
\alias{quant.trans}


\title{Non-parametric bootstrap resampling function in package `multtest'}

\description{Given a data set and a closure, which consists of a function for computing the test statistic and its enclosing environment, this function produces a non-parametric bootstrap estimated test statistics null distribution. The observations in the data are resampled using the ordinary non-parametric bootstrap is used to produce an estimated test statistics distribution. This distribution is then transformed to produce the null distribution. Options for transforming the nonparametric bootstrap distribution include \code{center.only}, \code{center.scale}, and \code{quant.trans}. Details are given below.  These functions are called by \code{MTP} and \code{EBMTP}.
}

\usage{
boot.null(X, label, stat.closure, W = NULL, B = 1000, test, nulldist, theta0 = 0, tau0 = 1, marg.null = NULL, marg.par = NULL, 
    ncp = 0, perm.mat, alternative = "two.sided", seed = NULL, 
    cluster = 1, dispatch = 0.05, keep.nulldist, keep.rawdist)

boot.resample(X, label, p, n, stat.closure, W, B, test)

center.only(muboot, theta0, alternative)

center.scale(muboot, theta0, tau0, alternative)

quant.trans(muboot, marg.null, marg.par, ncp, alternative, perm.mat) 

}

\arguments{
  \item{X}{A matrix, data.frame or ExpressionSet containing the raw data. In the case of an ExpressionSet, \code{exprs(X)} is the data of interest and \code{pData(X)} may contain outcomes and covariates of interest. For \code{boot.resample} \code{X} must be a matrix. For currently implemented tests, one hypothesis is tested for each row of the data.}
  \item{label}{A vector containing the class labels for t- and F-tests.}
  \item{stat.closure}{A closure for test statistic computation, like those produced internally by the \code{MTP} function. The closure consists of a function for computing the test statistic and its enclosing environment, with bindings for relevant additional arguments (such as null values, outcomes, and covariates).}
  \item{W}{A vector or matrix containing non-negative weights to be used in computing the test statistics. If a matrix, \code{W} must be the same dimension as \code{X} with one weight for each value in \code{X}. If a vector, \code{W} may contain one weight for each observation (i.e. column) of \code{X} or one weight for each variable (i.e. row) of \code{X}. In either case, the weights are duplicated appropriately. Weighted F-tests are not available. Default is 'NULL'.}
  \item{B}{The number of bootstrap iterations (i.e. how many resampled data sets) or the number of permutations (if \code{nulldist} is 'perm'). Can be reduced to increase the speed of computation, at a cost to precision. Default is 1000.}
  \item{test}{Character string specifying the test statistics to use. See \code{MTP} for a list of tests.}
  \item{theta0}{The value used to center the test statistics. For tests based on a form of t-statistics, this should be zero (default). For F-tests, this should be 1.}
  \item{tau0}{The value used to scale the test statistics. For tests based on a form of t-statistics, this should be 1 (default). For F-tests, this should be 2/(K-1), where K is the number of groups. This argument is missing when \code{center.only} is chosen for transforming the raw bootstrap test statistics.}
  \item{marg.null}{If \code{nulldist='boot.qt'}, the marginal null distribution to use for quantile transformation.  Can be one of 'normal', 't', 'f' or 'perm'.  Default is 'NULL', in which case the marginal null distribution is selected based on choice of test statistics.  Defaults explained below. If 'perm', the user must supply a vector or matrix of test statistics corresponding to another marginal null distribution, perhaps one created externally by the user, and possibly referring to empirically derived \emph{marginal permutation distributions}, although the statistics could represent any suitable choice of marginal null distribution.}   
  \item{marg.par}{If \code{nulldist='boot.qt'}, the parameters defining the marginal null distribution in \code{marg.null} to be used for quantile transformation.  Default is 'NULL', in which case the values are selected based on choice of test statistics and other available parameters (e.g., sample size, number of groups, etc.).  Defaults explained below.  User can override defaults, in which case a matrix of marginal null distribution parameters can be accepted.  Providing a matrix of values allows the user to perform multiple testing using parameters which may vary with each hypothesis, as may be desired in common-quantile minP procedures. In this way, factors affecting multiple testing procedure performance such as sample size or missingness may be assessed.}
  \item{ncp}{If \code{nulldist='boot.qt'}, a value for a possible noncentrality parameter to be used during marginal quantile transformation. Default is 'NULL'.}
  \item{perm.mat}{If \code{nulldist='boot.qt'} and \code{marg.null='perm'}, a  matrix of user-supplied test statistics from a particular distribution to be used during marginal quantile transformation.  The statistics may represent empirically derived marginal permutation values, may be theoretical values, or may represent a sample from some other suitable choice of marginal null distribution.}
  \item{alternative}{Character string indicating the alternative hypotheses, by default 'two.sided'. For one-sided tests, use 'less' or 'greater' for null hypotheses of 'greater than or equal' (i.e. alternative is 'less') and 'less than or equal', respectively.}
  \item{seed}{Integer or vector of integers to be used as argument to \code{set.seed} to set the seed for the random number generator for bootstrap resampling. This argument can be used to repeat exactly a test performed with a given seed. If the seed is specified via this argument, the same seed will be returned in the seed slot of the MTP object created. Else a random seed(s) will be generated, used and returned. Vector of integers used to specify seeds for each node in a cluster used to to generate a bootstrap null distribution.}
  \item{cluster}{Integer of 1 or a cluster object created through the package snow. With cluster=1, bootstrap is implemented on single node. Supplying a cluster object results in the bootstrap being implemented in parallel on the provided nodes. This option is only available for the bootstrap procedure.}
  \item{csnull}{DEPRECATED as of \code{multtest} v. 2.0.0 given expanded null distribution options. Previously, this argument was an indicator of whether the bootstrap estimated test statistics distribution should be centered and scaled (to produce a null distribution) or not. If \code{csnull=FALSE}, the (raw) non-null bootstrap estimated test statistics distribution was returned.  If the non-null bootstrap distribution should be returned, this object is now stored in the 'rawdist' slot when \code{keep.rawdist=TRUE}.} 
  \item{dispatch}{The number or percentage of bootstrap iterations to dispatch at a time to each node of the cluster if a computer cluster is used. If dispatch is a percentage, \code{B*dispatch} must be an integer. If dispatch is an integer, then \code{B/dispatch} must be an integer. Default is 5 percent.}
  \item{p}{An integer of the number of variables of interest to be tested.}
  \item{n}{An integer of the total number of samples.}
  \item{muboot}{A matrix of bootstrapped test statistics.}
  \item{keep.nulldist}{Logical indicating whether to return the computed bootstrap null distribution, by default 'TRUE'.  Not available for \code{nulldist}='perm'. Note that this matrix can be quite large.}
  \item{keep.rawdist}{Logical indicating whether to return the computed non-null (raw) bootstrap distribution, by default 'FALSE'.  Not available for when using \code{nulldist}='perm' or 'ic'.  Note that this matrix can become quite large.  If one wishes to use subsequent calls to \code{update} in which one updates choice of bootstrap null distribution, \code{keep.rawdist} must be TRUE.  To save on memory, \code{update} only requires that one of \code{keep.nulldist} or \code{keep.rawdist} be 'TRUE'.}
}

\value{
A list with the following elements:
  \item{rawboot}{If \code{keep.rawdist=TRUE}, the matrix of non-null, non-transformed bootstrap test statistics.  If 'FALSE', an empty matrix with dimension 0-by-0.}
  \item{muboot}{If \code{keep.rawdist=TRUE} (default), the matrix of appropriately transformed null test statistics as given by one of \code{center.scale}, \code{center.only}, or \code{quant.trans}. This is the estimated joint test statistics null distribution. \cr

Both list elements \code{rawboot} and \code{muboot} contain matrices of dimension the number of hypotheses (typically \code{nrow(X)}) by the number of bootstrap iterations (\code{B}). Each row of \code{muboot} is the bootstrap estimated marginal null distribution for a single hypothesis.  For \code{boot.null} and \code{center.scale}, each column of \code{muboot} is a centered and scaled resampled vector of test statistics.  For \code{boot.null} and \code{center.only}, each column of \code{muboot} is a centered, resampled vector of test statistics.\cr

For \code{boot.null} and \code{quant.trans}, each column of \code{muboot} is a marginal null quantile-transformed resampled vector of test statistics.  For each choice of marginal null distribution (defined by \code{marg.null} and \code{marg.par}), a random sample of size B is drawn and then rearranged based on the ranks of the marginal test statistics bootstrap distribution corresponding to each hypothesis (typically within rows of \code{X}).  This means that using \code{quant.trans} will set the RNG seed ahead by B * the number of hypotheses (similarly, typically \code{nrow(X)}).  Tie breaks in the marginal non-null bootstrap distribution are implemented inside the internal function \code{marg.samp} called by \code{quant.trans}.  Default values of \code{marg.null} and \code{marg.par} are available based on choice of test statistics, sample size 'n', and various other parameters.  By the time \code{boot.null} is called in either the \code{MTP} or \code{EBMTP} functions, the default marginal null distribution settings have already been formatted and passed in their correct form to \code{boot.null}.  These default values correspond to:

\describe{
\item{t.onesamp:}{t-distribution with df=n-1;}
\item{t.twosamp.equalvar:}{t-distribution with df=n-2;}
\item{t.twosamp.unequalvar:}{N(0,1);}
\item{t.pair:}{t-distribution with df=n-1, where n is the number of unique samples, i.e., the number of observed differences/paired samples;}
\item{f:}{F-distribution with df1=k-1, df2=n-k, for k groups;}
\item{f.block:}{NA. Only available with permutation distribution;}
\item{f.twoway:}{F-distribution with df1=k-1,df2=n-k*l, for k groups and l blocks;}
\item{lm.XvsZ:}{N(0,1);}
\item{lm.YvsXZ:}{N(0,1);}
\item{coxph.YvsXZ:}{N(0,1);}
\item{t.cor}{t-distribution with df=n-2;}
\item{z.cor}{N(0,1).}
}
The above defaults, however, can be overridden by manually setting values of \code{marg.null} and \code{marg.par}. \cr 

The \code{rawboot} and \code{muboot} objects are returned in the slots \code{rawdist} and \code{nulldist} of an object of class \code{MTP} or \code{EBMTP} when the arguments \code{keep.rawdist} or \code{keep.nulldist} to the \code{MTP} function are TRUE. For \code{boot.resample} a matrix of bootstrap samples prior to null transformation is returned.
}
}

\references{
M.J. van der Laan, S. Dudoit, K.S. Pollard (2004), Augmentation Procedures for Control of the Generalized Family-Wise Error Rate and Tail Probabilities for the Proportion of False Positives, Statistical Applications in Genetics and Molecular Biology, 3(1). 
\url{http://www.bepress.com/sagmb/vol3/iss1/art15/}

M.J. van der Laan, S. Dudoit, K.S. Pollard (2004), Multiple Testing. Part II. Step-Down Procedures for Control of the Family-Wise Error Rate, Statistical Applications in Genetics and Molecular Biology, 3(1).
\url{http://www.bepress.com/sagmb/vol3/iss1/art14/}

S. Dudoit, M.J. van der Laan, K.S. Pollard (2004), Multiple Testing. Part I. Single-Step Procedures for Control of General Type I Error Rates, Statistical Applications in Genetics and Molecular Biology, 3(1).
\url{http://www.bepress.com/sagmb/vol3/iss1/art13/}

Katherine S. Pollard and Mark J. van der Laan, "Resampling-based Multiple Testing: Asymptotic Control of Type I Error and Applications to Gene Expression Data" (June 24, 2003). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 121.
\url{http://www.bepress.com/ucbbiostat/paper121}

M.J. van der Laan and A.E. Hubbard (2006), Quantile-function Based Null Distributions in Resampling Based Multiple Testing, Statistical Applications in Genetics and Molecular Biology, 5(1).
\url{http://www.bepress.com/sagmb/vol5/iss1/art14/}

S. Dudoit and M.J. van der Laan.  Multiple Testing Procedures and Applications to Genomics.  Springer Series in Statistics. Springer, New York, 2008.

}

\author{Katherine S. Pollard, Houston N. Gilbert, and Sandra Taylor, with design contributions from Sandrine Dudoit and Mark J. van der Laan.}

\note{Thank you to Duncan Temple Lang and Peter Dimitrov for suggestions about the code.}

\seealso{\code{\link{corr.null}}, \code{\link{MTP}}, \code{\link{MTP-class}}, \code{\link{EBMTP}}, \code{\link{EBMTP-class}}, \code{\link{get.Tn}}, \code{\link{ss.maxT}}, \code{\link{mt.sample.teststat}},\code{\link{get.Tn}}, \code{\link{wapply}}, \code{\link{boot.resample}}}

\examples{

set.seed(99)
data<-matrix(rnorm(90),nr=9)

#closure
ttest<-meanX(psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE)

#test statistics
obs<-get.Tn(X=data,stat.closure=ttest,W=NULL)

#bootstrap null distribution (B=100 for speed, default nulldist, "boot.cs")
nulldistn<-boot.null(X=data,W=NULL,stat.closure=ttest,B=100,test="t.onesamp",
	nulldist="boot.cs",theta0=0,tau0=1,alternative="two.sided",
	keep.nulldist=TRUE,keep.rawdist=FALSE)$muboot

#bootstrap null distribution with marginal quantile transformation showing
#default values that are passed to marg.null and marg.par arguments
nulldistn.qt<-boot.null(X=data,W=NULL,stat.closure=ttest,B=100,test="t.onesamp",
	nulldist="boot.qt",theta0=0,tau0=1,alternative="two.sided",
	keep.nulldist=TRUE,keep.rawdist=FALSE,marg.null="t",
	marg.par=matrix(9,nr=10,nc=1))$muboot

#unadjusted p-values
rawp<-apply((obs[1,]/obs[2,])<=nulldistn,1,mean)
sum(rawp<=0.01)

rawp.qt<-apply((obs[1,]/obs[2,])<=nulldistn.qt,1,mean)
sum(rawp.qt<=0.01)
}

\keyword{manip}
\keyword{internal}



