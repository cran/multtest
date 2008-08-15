\name{meanX}
\alias{meanX}
\alias{diffmeanX}
\alias{FX}
\alias{blockFX}
\alias{twowayFX}
\alias{lmX}
\alias{lmY}
\alias{coxY}
\alias{get.Tn}

\title{Functions to create test statistic closures and apply them to data}

\description{
The package \code{multtest} uses closures in the function \code{MTP} to compute test statistics. The closure used depends on the value of the argument \code{test}. These functions create the closures for different tests, given any additional variables, such as outcomes or covariates. The function \code{get.Tn} calls \code{wapply} to apply one of these closures to observed data (and possibly weights).
}

\usage{
meanX(psi0 = 0, na.rm = TRUE, standardize = TRUE, 
alternative = "two.sided", robust = FALSE)

diffmeanX(label, psi0 = 0, var.equal = FALSE, na.rm = TRUE, 
standardize = TRUE, alternative = "two.sided", robust = FALSE)

FX(label, na.rm = TRUE, robust = FALSE)

blockFX(label, na.rm = TRUE, robust = FALSE)

twowayFX(label, na.rm = TRUE, robust = FALSE)

lmX(Z = NULL, n, psi0 = 0, na.rm = TRUE, standardize = TRUE, 
alternative = "two.sided", robust = FALSE)

lmY(Y, Z = NULL, n, psi0 = 0, na.rm = TRUE, standardize = TRUE, 
alternative = "two.sided", robust = FALSE)

coxY(surv.obj, strata = NULL, psi0 = 0, na.rm = TRUE, standardize = TRUE, 
alternative = "two.sided", init = NULL, method = "efron")

get.Tn(X, stat.closure, W = NULL)
}

\arguments{
  \item{X}{A matrix, data.frame or ExpressionSet containing the raw data. In the case of an ExpressionSet, \code{exprs(X)} is the data of interest and \code{pData(X)} may contain outcomes and covariates of interest. For currently implemented tests, one hypothesis is tested for each row of the data.}
  \item{W}{A vector or matrix containing non-negative weights to be used in computing the test statistics. If a matrix, \code{W} must be the same dimension as \code{X} with one weight for each value in \code{X}. If a vector, \code{W} may contain one weight for each observation (i.e. column) of \code{X} or one weight for each variable (i.e. row) of \code{X}. In either case, the weights are duplicated apporpraiately. Weighted f-tests are not available. Default is 'NULL'.}
  \item{label}{A vector containing the class labels for t- and f-tests. For the \code{blockFX} function, observations are divided into \code{l} blocks of \code{n/l} observations. Within each block there may be \code{k} groups with \code{k>2}. For this test, there is only one observation per block*group combination. The labels (and corresponding rows of \code{Z} and columns of \code{X} and \code{W}) must be ordered by block and within each block ordered by group. Groups must be labeled with integers \code{1,...,k}. For the \code{twowayFX} function, observations are divided into \code{l} blocks. Within each block there may be \code{k} groups with \code{k>2}. There must be more than one observation per group*block combination for this test. The labels (and corresponding rows of \code{Z} and columns of \code{X} and \code{W}) must be ordered by block and within each block ordered by group. Groups must be labeled with integers \code{1,...,k}.}
  \item{Y}{A vector or factor containing the outcome of interest for linear models. This may be a continuous or polycotomous dependent variable.}
  \item{surv.object}{A survival object as returned by the \code{Surv} function, to be used as response in \code{coxY}.} 
  \item{Z}{A vector, factor, or matrix containing covariate data to be used in the linear regression models. Each variable should be in one column.}
  \item{strata}{A vector, factor, or matrix containing covariate data to be used in the Cox regression models. Covariate data will be converted to a factor variable (via the \code{strata} function) for use in the \code{coxph} function. Each variable should be in one column.} 
  \item{n}{The sample size, e.g. \code{length(Y)} or \code{nrow(Z)}.}
  \item{psi0}{Hypothesized null value for the parameter of interest (e.g. mean or difference in means), typically zero (default).}
  \item{var.equal}{Indicator of whether to use t-statistics that assume equal variance in the two groups when computing the denominator of the test statistics.}
  \item{na.rm}{Logical indicating whether to remove observations with an NA. Default is 'TRUE'.}
  \item{standardize}{Logical indicating whether to use the standardized version of the test statistics (usual t-statistics are standardized). Default is 'TRUE'.}
  \item{alternative}{Character string indicating the alternative hypotheses, by default 'two.sided'. For one-sided tests, use 'less' or 'greater' for null hypotheses of 'greater than or equal' (i.e. alternative is 'less') and 'less than or equal', respectively.}
  \item{robust}{Logical indicating whether to use robust versions of the test statistics.}
  \item{init}{Vector of initial values of the iteration in \code{coxY} function, as used in \code{coxph} in the \code{survival} package. Default initial value is zero for all variables (\code{init=NULL}).}
  \item{method}{A character string specifying the method for tie handling in \code{coxY} function, as used in \code{coxph} in the \code{survival} package. Default is "efron".}

}

\details{
The use of closures, in the style of the \code{genefilter} package, allows uniform data input for all MTPs and facilitates the extension of the package's functionality by adding, for example, new types of test statistics. 
Specifically, for each value of the \code{MTP} argument \code{test}, a closure is defined which consists of a function for computing the test statistic (with only two arguments, a data vector \code{x} and a corresponding weight vector \code{w}, with default value of \code{NULL}) and its enclosing environment, with bindings for relevant additional arguments. These arguments may include null values \code{psi0}, outcomes (\code{Y}, \code{label}, \code{surv.object}), and covariates \code{Z}. The vectors \code{x} and \code{w} are rows of the matrices \code{X} and \code{W}.

In the \code{MTP} function, the closure is first used to compute the vector of observed test statistics, and then, in each bootstrap iteration, to produce the estimated joint null distribution of the test statistics. In both cases, the function \code{get.Tn} is used to apply the closure to rows of the matrices of data (\code{X}) and weights (\code{W}). Thus, new test statistics can be added to \code{multtest} package by simply defining a new closure and adding a corresponding value for the \code{test} argument to the \code{MTP} function.
}

\value{
For \code{meanX}, \code{diffmeanX}, \code{FX}, \code{blockFX}, \code{twowayFX}, \code{lmX}, \code{lmY}, and \code{coxY}, a closure consisting of a function for computing test statistics and its enclosing environment. For \code{get.Tn}, the observed test statistics stored in a matrix \code{obs} with numerator (possibly absolute value or negative, depending on the value of alternative) in the first row, denominator in the second row, and a 1 or -1 in the third row (depending on the value of alternative). The vector of observed test statistics is obs[1,]*obs[3,]/obs[2,].
}

\author{Katherine S. Pollard, Sandra Taylor, with design contributions from Duncan Temple Lang, Sandrine Dudoit and Mark J. van der Laan}

\seealso{\code{\link{MTP}}, \code{\link{get.Tn}}, \code{\link{wapply}}, \code{\link{boot.resample}}}

\examples{
data<-matrix(rnorm(200),nr=20)
ttest<-meanX(psi0=0,na.rm=TRUE,standardize=TRUE,alternative="two.sided",robust=FALSE)
obs<-wapply(data,1,ttest,W=NULL)
statistics<-obs[1,]*obs[3,]/obs[2,]
statistics


FData <- matrix(rnorm(5*60),nr=5)
label<-rep(c(rep(1,10), rep(2,10), rep(3,10)),2)
twowayf<-twowayFX(label)
obs<-wapply(FData,1,twowayf,W=NULL)
statistics<-obs[1,]*obs[3,]/obs[2,]
statistics
}

\keyword{htest}
\keyword{internal}



