\name{get.index}

\alias{get.index}

\title{Function to compute indices for ordering hypotheses in Package 'multtest'}

\description{
The hypotheses tested in a multiple testing procedure (MTP), can be ordered based on the output of that procedure. This function orders hypotheses based on adjusted p-values, then unadjusted p-values (to break ties in adjusted p-values), and finally test statistics (to break remaining ties).
}

\usage{
get.index(adjp, rawp, stat)
}

\arguments{
  \item{adjp}{Numeric vector of adjusted p-values.}
  \item{rawp}{Numeric vector of unadjusted ("raw") marginal p-values.}
  \item{stat}{Numeric vector of test statistics.}
}

\value{
Numeric vector of indices so that the hypotheses can be ordered accroding to significance (smallest p-values and largest test statistics first). This function is used in the plot method for objects of class \code{MTP} to order adjusted p-values for graphical summaries. The summary method for objects of class \code{MTP} will return these indices as its second component.
}

\author{Katherine S. Pollard, \url{http://lowelab.ucsc.edu/katie/}}

\seealso{\code{\link{MTP}}, \code{\link{plot,MTP-method}}, \code{\link{summary,MTP-method}}}

\examples{
data<-matrix(rnorm(200),nr=20)
mtp<-MTP(X=data,test="t.onesamp")
index<-get.index(adjp=mtp@adjp,rawp=mtp@rawp,stat=mtp@statistic)
mtp@statistic[index]
mtp@estimate[index]
apply(data[index,],1,mean)
}

\keyword{htest}
\keyword{internal}






