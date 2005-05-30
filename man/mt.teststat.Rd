\name{mt.teststat}
\alias{mt.teststat}
\alias{mt.teststat.num.denum}

\title{Computing test statistics for each row of a data frame}
\usage{
mt.teststat(X,classlabel,test="t",na=.mt.naNUM,nonpara="n")
mt.teststat.num.denum(X,classlabel,test="t",na=.mt.naNUM,nonpara="n")
}
\description{
  These functions provide a convenient way to compute test statistics,
  e.g., two-sample Welch t-statistics, Wilcoxon statistics,
  F-statistics, paired t-statistics,
  block F-statistics, for each row of a data frame. 
}
\arguments{
 \item{X}{A data frame or matrix, with \eqn{m} rows corresponding to variables
    (hypotheses) and\eqn{n} columns to observations. In the case of gene 
    expression data, rows
    correspond to genes and columns to mRNA samples. The data can
    be read using \code{\link{read.table}}.
  }
  \item{classlabel}{
     A vector of integers corresponding to observation (column)
    class labels. For \eqn{k} classes, the labels must be integers
    between 0 and \eqn{k-1}. For the \code{blockf} test option,
    observations may be divided into
    \eqn{n/k} blocks of \eqn{k} observations each. The observations are
    ordered by block, and within each block, they are labeled using the
    integers 0 to \eqn{k-1}.
  }	
  \item{test}{A character string specifying the statistic to be
    used to test the null hypothesis of no association between the
    variables and the class labels.\cr
    If \code{test="t"}, the tests are based on two-sample Welch t-statistics
    (unequal variances).  \cr
    If \code{test="t.equalvar"}, the tests are based on two-sample
    t-statistics with equal variance for the two samples. The
    square of the t-statistic is equal to an F-statistic for \eqn{k=2}. \cr
    If \code{test="wilcoxon"}, the tests are based on standardized rank sum Wilcoxon statistics.\cr
    If \code{test="f"}, the tests are based on F-statistics.\cr
    If \code{test="pairt"}, the tests are based on paired t-statistics. The
    square of the paired t-statistic is equal to a block F-statistic for \eqn{k=2}. \cr
    If \code{test="blockf"}, the tests are based on F-statistics which
    adjust for block differences
    (cf. two-way analysis of variance).
  }
  \item{na}{Code for missing values (the default is \code{.mt.naNUM=--93074815.62}).
    Entries with missing values will be ignored in the computation,
    i.e., test statistics will be based on a smaller sample size. This
    feature has not yet fully implemented.
  }
  \item{nonpara}{If \code{nonpara}="y", nonparametric test statistics are computed based on ranked data. \cr
    If  \code{nonpara}="n", the original data are used.}
}

\value{
  For \code{\link{mt.teststat}}, a vector of test statistics for each row (gene). \cr \cr
  For \code{\link{mt.teststat.num.denum}}, a data frame with \cr
  \item{teststat.num}{the numerator of the test statistics for each row, depending on the
    specific \code{test} option.}
  \item{teststat.denum}{the denominator of the test statistics for each row, depending on the
    specific \code{test} option.}
  }


\author{Yongchao Ge, \email{yongchao.ge@mssm.edu}, \cr
Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{mt.maxT}}, \code{\link{mt.minP}}, \code{\link{golub}}.}

\examples{
# Gene expression data from Golub et al. (1999)
data(golub)

teststat<-mt.teststat(golub,golub.cl)
qqnorm(teststat)
qqline(teststat)

tmp<-mt.teststat.num.denum(golub,golub.cl,test="t")
num<-tmp$teststat.num
denum<-tmp$teststat.denum
plot(sqrt(denum),num)

tmp<-mt.teststat.num.denum(golub,golub.cl,test="f")

}
		
\keyword{univar}

	

	
