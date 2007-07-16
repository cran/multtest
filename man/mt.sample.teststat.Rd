\name{mt.sample.teststat}
\title{Permutation distribution of test statistics and raw (unadjusted) p-values}
\alias{mt.sample.teststat}
\alias{mt.sample.rawp}
\alias{mt.sample.label}
\usage{
mt.sample.teststat(V,classlabel,test="t",fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
mt.sample.rawp(V,classlabel,test="t",side="abs",fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
mt.sample.label(classlabel,test="t",fixed.seed.sampling="y",B=10000)
}
\description{
  These functions provide tools to investigate the permutation distribution
  of test statistics, raw (unadjusted) \eqn{p}-values, and class labels.
}
\arguments{
  \item{V}{A numeric vector containing the data for one of the variables (genes).}

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
  \item{side}{A character string specifying the type of rejection region.\cr
    If \code{side="abs"}, two-tailed tests, the null hypothesis is rejected for large absolute values of the test statistic.\cr
    If \code{side="upper"}, one-tailed tests, the null hypothesis is rejected for large values of the test statistic.\cr
    If \code{side="lower"}, one-tailed tests,  the null hypothesis is rejected for small values of the test statistic.
  }
  \item{fixed.seed.sampling}{If \code{fixed.seed.sampling="y"}, a
    fixed seed sampling procedure is used, which may double the
    computing time, but will not use extra memory to store the
    permutations. If \code{fixed.seed.sampling="n"}, permutations will
    be stored in memory.  For the \code{blockf} test, the option \code{n} was not implemented as it requires too much memory.
  }
  \item{B}{The number of permutations. For a complete
    enumeration, \code{B} should be 0 (zero) or any number not less than
    the total number of permutations.
  }
  \item{na}{Code for missing values (the default is \code{.mt.naNUM=--93074815.62}).
    Entries with missing values will be ignored in the computation,
    i.e., test statistics will be based on a smaller sample size. This
    feature has not yet fully implemented.
  }
  \item{nonpara}{If \code{nonpara}="y", nonparametric test statistics are computed based on ranked data. \cr
    If  \code{nonpara}="n", the original data are used.
  }

}
\value{
  For \code{\link{mt.sample.teststat}},  a vector containing \code{B} permutation test statistics. \cr \cr
  For \code{\link{mt.sample.rawp}},  a vector containing \code{B} permutation unadjusted \eqn{p}-values. \cr\cr 
  For \code{\link{mt.sample.label}}, a matrix containing \code{B}
  sets of permuted class labels. Each row corresponds to one permutation.
}

\examples{

# Gene expression data from Golub et al. (1999)
data(golub)

mt.sample.label(golub.cl,B=10)

permt<-mt.sample.teststat(golub[1,],golub.cl,B=1000)
qqnorm(permt)
qqline(permt)

permt<-mt.sample.teststat(golub[50,],golub.cl,B=1000)
qqnorm(permt)
qqline(permt)

permp<-mt.sample.rawp(golub[1,],golub.cl,B=1000)
hist(permp)
}

\author{Yongchao Ge, \email{yongchao.ge@mssm.edu}, \cr
Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{mt.maxT}}, \code{\link{mt.minP}}, \code{\link{golub}}.}

\keyword{manip}

