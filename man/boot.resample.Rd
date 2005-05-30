\name{boot.resample}
\alias{boot.resample}

\title{Non-parametric bootstrap resampling function in package `multtest'}

\description{Given a data set and a closure, which consists of a function for computing the test statistic and its enclosing environment, this function produces a non-parametric bootstrap estimated test statistics null distribution. The observations in the data are resampled using the ordinary non-parametric bootstrap and used to produce an estimated test statistics distribution. This distribution is then centered and scaled to produce the null distribution. This function is called by \code{MTP}.
}

\usage{
boot.resample(X, stat.closure, W = NULL, B = 1000, theta0 = 0, tau0 = 1)
}

\arguments{
  \item{X}{A matrix or data.frame containing the raw data. For currently implemented tests, one hypothesis is tested for each row of the data.}
  \item{stat.closure}{A closure for test statistic computation, like those produced internally by the \code{MTP} function. The closure consists of a function for computing the test statistic and its enclosing environment, with bindings for relevant additional arguments (such as null values, outcomes, and covariates).}
  \item{W}{A vector or matrix containing non-negative weights to be used in computing the test statistics. If a matrix, \code{W} must be the same dimension as \code{X} with one weight for each value in \code{X}. If a vector, \code{W} may contain one weight for each observation (i.e. column) of \code{X} or one weight for each variable (i.e. row) of \code{X}. In either case, the weights are duplicated appropriately. Weighted f-tests are not available. Default is 'NULL'.}
  \item{B}{The number of bootstrap iterations (i.e. how many resampled data sets) or the number of permutations (if \code{nulldist} is 'perm'). Can be reduced to increase the speed of computation, at a cost to precision. Default is 1000.}
  \item{theta0}{The value used to center the test statistics. For tests based on a form of t-statistics, this should be zero (default). For f-tests, this should be 1.}
  \item{tau0}{The value used to scale the test statistics. For tests based on a form of t-statistics, this should be 1 (default). For f-tests, this should be 2/(K-1), where K is the number of groups.}
}

\value{
A matrix of dimension number of hypotheses (\code{nrow(X)}) by number of bootstrap iterations (\code{B}). This is the estimated joint test statistics null distribution. Each column is a centered and scaled resampled vector of test statistics. Each row is the bootstrap estimated marginal null distribution for a single hypothesis. This object is returned in slot \code{nulldist} of an object of class \code{MTP} when the argument \code{keep.nulldist} to the \code{MTP} function is TRUE.
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
}

\author{Katherine S. Pollard, \url{http://lowelab.ucsc.edu/katie/}\cr
with design contributions from Sandrine Dudoit and Mark J. van der Laan.}

\notes{Thank you to Peter Dimitrov for suggestions about the code.}

\seealso{\code{\link{MTP}}, \code{\link{MTP-class}}, \code{\link{get.Tn}}, \code{\link{ss.maxT}}, \code{\link{mt.sample.teststat}}}

\examples{

#data example: ALL data set
set.seed(99)
data<-matrix(rnorm(90),nr=9)

#closure
ttest<-meanX(psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE)

#test statistics
obs<-get.Tn(X=data,stat.closure=ttest,W=NULL)

#bootstrap null distribution (B=100 for speed)
nulldistn<-boot.resample(X=data,W=NULL,stat.closure=ttest,B=100,theta0=0,tau0=1)

#unadjusted p-values
rawp<-apply((obs[1,]/obs[2,])<=nulldistn,1,mean)
sum(rawp<=0.01)

}

\keyword{manip}
\keyword{internal}



