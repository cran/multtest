\name{corr.null}
\alias{corr.null}
\alias{tQuantTrans}

\title{Function to estimate a test statistics joint null distribution for t-statistics via the vector influence curve}

\description{For a broad class of testing problems, such as the test of single-parameter null hypotheses using t-statistics, a proper, asymptotically valid test statistics joint null distribution is the multivariate Gaussian distribution with mean vector zero and covariance matrix equal to the correlation matrix of the vector influence curve for the estimator of the parameter of interest.  The function \code{corr.null} estimates the correlation matrix of the vector influence curve for such parameters and returns samples from the corresponding normal distribution.  Arguments to the function allow for refinements in calculating the resulting null distribution estimate.} 


\usage{
corr.null(X, W = NULL, Y = NULL, Z = NULL, test = "t.twosamp.unequalvar", 
    alternative = "two-sided", use = "pairwise", B = 1000, MVN.method = "mvrnorm", 
    penalty = 1e-06, ic.quant.trans = FALSE, marg.null = NULL, 
    marg.par = NULL, perm.mat = NULL) 
}

\arguments{
  \item{X}{A matrix, data.frame or ExpressionSet containing the raw data. In the case of an ExpressionSet, \code{exprs(X)} is the data of interest and \code{pData(X)} may contain outcomes and covariates of interest. For most currently implemented tests (exception: tests involving correlation parameters), one hypothesis is tested for each row of the data.}
  \item{W}{A matrix containing non-negative weights to be used in computing the test statistics.  Must be same dimension as \code{X}.}
  \item{Y}{A vector, factor, or \code{Surv} object containing the outcome of interest.}
  \item{Z}{A vector, factor, or matrix containing covariate data to be used in linear regression models. Each variable should be in one column, so that \code{nrow(Z)=ncol(X)}. By the time the function is called, this argument contains a 'design matrix' with the variable to be tested in the first column, additional covariates in the remaining columns, and no intercept column.}
  \item{test}{Character string specifying the test statistics to use, by default 't.twosamp.unequalvar'. See details (below) for a list of tests.}
  \item{alternative}{Character string indicating the alternative hypotheses, by default 'two.sided'. For one-sided tests, use 'less' or 'greater' for null hypotheses of 'greater than or equal' (i.e. alternative is 'less') and 'less than or equal', respectively.}
  \item{use}{Similar to the options in \code{cor}, a character string giving a method for computing covariances in the presence of missing values.  Default is 'pairwise', which allows for the covariance/correlation matrix to be calculated using the most information possible when \code{NA}s are present.} 
  \item{B}{The number of samples to be drawn from the normal distribution. Default is 1000.}
  \item{MVN.method}{Character string of either of 'mvrnorm' or 'Cholesky' designating how correlated normal test statistics are to be generated.  Selecting 'mvrnorm' uses the function of the same name found in the \code{MASS} library, whereas 'Cholesky' relies on a Cholesky decomposition. Default is 'mvrnorm'.}
  \item{penalty}{If \code{MVN.method='Cholesky'}, the value in \code{penalty} is added to all diagonal elements of the estimated test statistics correlation matrix to ensure that the matrix is positive definite and that internal calls to \code{'chol'} do not return an error.  Default is 1e-6.}
  \item{ic.quant.trans}{A logical indicating whether or not a marginal quantile transformation using a t-distribution or user-supplied marginal distribution (stored in \code{perm.mat}) should be applied to the multivariate normal null distribution.  Defaults for \code{marg.null} and \code{marg.par} exist, but can also be specified by the user (see below). Default is 'FALSE'.}
  \item{marg.null}{If \code{ic.quant.trans=TRUE}, a character string naming the marginal null distribution to use for quantile transformation.  Can be one of, 't' or 'perm'.  Default is 'NULL', in which case the marginal null distribution is selected based on choice of test statistics.  Defaults explained below. If 'perm', the user must supply a vector or matrix of test statistics corresponding to another marginal null distribution, perhaps one created externally by the user, and possibly referring to empirically derived \emph{marginal permutation distributions}, although the statistics could represent any suitable choice of marginal null distribution.}   
  \item{marg.par}{If \code{ic.quant.trans=TRUE}, the parameters defining the marginal null distribution in \code{marg.null} to be used for quantile transformation.  Default is 'NULL', in which case the values are selected based on choice of test statistics and other available parameters (e.g., sample size, number of groups, etc.).  Defaults explained below.  User can override defaults, in which case a matrix of marginal null distribution parameters must be provided.  Providing a matrix allows the user to perform multiple testing using parameters which may vary with each hypothesis, as may be desired in common-quantile minP procedures}
  \item{perm.mat}{If \code{ic.quant.trans=TRUE},  a matrix of user-supplied test statistics from a particular distribution to be used during marginal quantile transformation.  Supplying a vector of test statistics will apply the same vector to each hypothesis.  The statistics may represent empirically derived marginal permutation values, may be theoretical values, or may represent a sample from some other suitable choice of marginal null distribution.}
}

\details{
This function is called internally when the argument \code{nulldist='ic'} is evaluated in the main user-level functions \code{MTP} or \code{EBMTP}.  Formatting of the data objects \code{X}, \code{W}, \code{Y}, and especially \code{Z} occurs at execution begin of the main user-level functions.\cr

Based on the value of \code{test}, the appropriate correlation matrix of the vector influence curve is calculated. Once the correlation matrix is obtained, one may sample vectors of null test statistics directly from a multivariate normal distribution rather than relying on permutation-based or bootstrap-based resampling.  Because the Gaussian distribution is continuous, we expect this choice of null distribution to suffer less from discreteness than either the permutation or the bootstrap distribution.  Additionally, in large-scale settings, use of null distributions derived from the vector influence function typically reduce computational bottlenecks associated with resampling methods.\cr

Because the influence curve null distributions have been implemented for parametric, standardized t-statistics, the options \code{robust} and \code{standardize} are not allowed. Influence curve null distributions are available for the following values of \code{test}: 't.onesamp', 't.pair', 't.twosamp.equalvar', 't.twosamp.unequalvar', 'lm.XvsZ', 'lm.YvsXZ', 't.cor', and 'z.cor'.\cr

In the simpler cases involving one-sample and two-sample tests of means, the correlation matrices are obtained via calls to \code{cor}.  For two-sample tests, the correlation matrix corresponds to the following transformation of the group-specific covariance matrices: cov(X(group1))/n1 + cov(X(group2))/n2, where n1 and n2 are sample sizes of each group. When weights are present, the internal function \code{IC.CorXW.NA} is called to calculate weighted estimates of the (group) covariance matrices from each subject's estimated vector influence curve.  The calculations are similar in spirit to those in \code{cov.wt}, but they are done in a way which allows for handling \code{NA} elements in the estimated vector influence curve IC_n.  The correlation matrix corresponding to IC_n * (IC_n)^t is calculated. \cr 

For linear regression models, \code{corr.null} calculates the vector influence curve associated
with each subject/sample.  The vector has length equal to the number of hypotheses.  The internal function \code{IC.Cor.NA} is used to calculate IC_n * (IC_n)^t in a manner which allows for NA-handling when the influence curve may contain missing elements.  For linear regression models of the form E[Y|X], IC_n takes the form (E[((X^t)X)^(-1)] (X^t)_i Y_i) - Y_i-hat.  Influence curves for correlation parameters are more complicated, and the user is referred to the references below.\cr

Once the correlation matrix sigma' corresponding to the variance covariance matrix of the vector influence curve sigma =IC_n * (IC_n)^t is obtained, one may sample from N(0,sigma') to obtain null test statistics.\cr

If \code{ic.quant.trans=TRUE}, the matrix of null test statistics can be quantile transformed to produce a matrix which accounts for the joint dependencies between test statistics (down columns), but which has marginal t-distributions (across rows).  If \code{marg.null} and \code{marg.par} are not specified (=NULL), the following default t-distributions are applied:\cr

\describe{
\item{t.onesamp}{df=n-1;}
\item{t.pair}{df=n-1, where n is the number of unique samples, i.e., the number of observed differences between paired samples;}
\item{t.twosamp.equalvar}{df=n-2;}
\item{t.twosamp.unequalvar}{df=n-1; N.B., this is not recommended, since the effective degrees of freedom are unknown.  With sufficiently large n, a normal approximation should yield similar results.}
\item{lm.XvsZ}{df=n-p, where p is the number of variables in the regression equation;}
\item{lm.YvsXZ}{df=n-p, where p is the number of variables in the regression equation;}
\item{t.cor}{df=n-2;}
\item{z.cor}{N.B., also not recommended.  Fisher's z-statistics are already normally distributed.  Marginal transformation to a t-distribution makes little sense.}
}
}

\value{
A matrix of null test statistics with dimension the number of hypotheses (typically \code{nrow(X)}) by the number of desired samples (\code{B}).
}

\references{
K.S. Pollard and Mark J. van der Laan, "Resampling-based Multiple Testing: Asymptotic Control of Type I Error and Applications to Gene Expression Data" (June 24, 2003). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 121.
\url{http://www.bepress.com/ucbbiostat/paper121}

S. Dudoit and M.J. van der Laan.  Multiple Testing Procedures and Applications to Genomics.  Springer Series in Statistics. Springer, New York, 2008.

H.N. Gilbert, M.J. van der Laan, and S. Dudoit, "Joint Multiple Testing Procedures for Inferring Genetic Networks from Lower-Order Conditional Independence Graphs" (2009). \emph{In preparation.}
}

\author{Houston N. Gilbert}

\seealso{\code{\link{boot.null}},\code{\link{MTP}}, \code{\link{MTP-class}}, \code{\link{EBMTP}}, \code{\link{EBMTP-class}}, \code{\link{get.Tn}}, \code{\link{ss.maxT}}, \code{\link{mt.sample.teststat}},\code{\link{get.Tn}}, \code{\link{wapply}}, \code{\link{boot.resample}}}

\examples{
set.seed(99)
data <- matrix(rnorm(10*50),nr=10,nc=50)
nulldistn.mvrnorm <- corr.null(data,t="t.onesamp",alternative="greater",B=5000)
nulldistn.chol <- corr.null(data,t="t.onesamp",MVN.method="Cholesky",penalty=1e-9)
nulldistn.t <- corr.null(data,t="t.onesamp",ic.quant.trans=TRUE)
dim(nulldistn.mvrnorm)

}

\keyword{htest}
\keyword{internal}
