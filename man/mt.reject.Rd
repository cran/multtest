\name{mt.reject}
\alias{mt.reject}
\title{Identity and number of rejected hypotheses }
\description{This function returns the identity and number of rejected hypotheses for several multiple testing procedures and different nominal Type I error rates.
}
\usage{
mt.reject(adjp, alpha)
}
\arguments{
  \item{adjp}{A matrix of adjusted \emph{p}-values, with rows
    corresponding to hypotheses and columns to multiple testing
    procedures. This matrix could be obtained from the function
    \code{\link{mt.rawp2adjp}}
    .}
  \item{alpha}{A vector of nominal Type I error rates.}
}
\value{
  A list with components
  \item{r}{A matrix containing the number of rejected hypotheses for several multiple testing procedures and different nominal Type I error rates. Rows correspond to Type I error rates and columns to multiple testing procedures.}
  \item{which}{A matrix of indicators for the rejection of individual hypotheses by different multiple testing procedures for a nominal Type I error rate \code{alpha[1]}. Rows correspond to hypotheses and columns to multiple testing procedures.}
}

\author{
  Sandrine Dudoit,  \url{http://www.stat.berkeley.edu/~sandrine}, \cr
  Yongchao Ge, \email{yongchao.ge@mssm.edu}.
}


\seealso{\code{\link{mt.maxT}}, \code{\link{mt.minP}}, \code{\link{mt.rawp2adjp}}, \code{\link{golub}}.}

\examples{
# Gene expression data from Golub et al. (1999)
# To reduce computation time and for illustrative purposes, we condider only
# the first 100 genes and use the default of B=10,000 permutations.
# In general, one would need a much larger number of permutations
# for microarray data.

data(golub)
smallgd<-golub[1:100,] 
classlabel<-golub.cl

# Permutation unadjusted p-values and adjusted p-values for maxT procedure
res<-mt.maxT(smallgd,classlabel)
mt.reject(cbind(res$rawp,res$adjp),seq(0,1,0.1))$r

}

\keyword{htest}
