\name{wapply}

\alias{wapply}

\title{Weighted version of the apply function}

\description{
A function to perform 'apply' on an matrix of data and corresponding matrix of weights. 
}

\usage{
wapply(X, MARGIN, FUN, W = NULL, ...)
}

\arguments{
  \item{X}{A matrix of data.}
  \item{MARGIN}{A vector giving the subscripts which the function will be
          applied over. 1 indicates rows, 2 indicates columns.}
  \item{FUN}{The function to be applied. In the case of functions like
          \code{+} the function name must be quoted.}
  \item{W}{An optional matrix of weights. When \code{W=NULL}, the usual \code{apply} function is called.}
  \item{\dots}{optional arguments to \code{FUN}.}
}

\details{
When weights are provided, these are passed to \code{FUN} along with the data \code{X}. For example, if \code{FUN=meanX}, each data value is multiplied by the corresponding weight before the mean is applied.
}

\value{
If each call to \code{FUN} returns a vector of length \code{n}, then \code{wapply} returns an array of dimension \code{c(n, dim(X)[MARGIN])} if \code{n > 1}. If \code{n = 1}, \code{wapply} returns a vector if \code{MARGIN} has length 1 and an array of dimension \code{dim(X)[MARGIN]} otherwise. If \code{n = 0}, the result has length 0 but not necessarily the "correct" dimension.

If the calls to \code{FUN} return vectors of different lengths, \code{wapply} returns a list of length \code{dim(X)[MARGIN]}.

This function is used in the package \code{multtest} to compute weighted versions of test statistics. It is called by the function \code{get.Tn} inside the user-level function \code{MTP}.
}

\author{Katherine S. Pollard}

\seealso{\code{\link{get.Tn}}, \code{\link{MTP}}}

\examples{
data<-matrix(rnorm(200),nr=20)
weights<-matrix(rexp(200,rate=0.1),nr=20)
wapply(X=data,MARGIN=1,FUN=mean,W=weights)
}

\keyword{internal}



