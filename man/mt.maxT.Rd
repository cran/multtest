\name{mt.maxT}
\alias{mt.maxT}
\alias{mt.minP}
\title{
  Step-down maxT and minP multiple testing procedures
}
\description{These functions compute permutation adjusted \eqn{p}-values for step-down multiple testing procedures described in Westfall & Young (1993).
}
\usage{
mt.maxT(X,classlabel,test="t",side="abs",fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
mt.minP(X,classlabel,test="t",side="abs",fixed.seed.sampling="y",B=10000,na=.mt.naNUM,nonpara="n")
}

\arguments{
  \item{X}{A data frame or matrix, with \eqn{m} rows corresponding to variables
    (hypotheses) and
    \eqn{n} columns to observations. In the case of gene expression data, rows
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

\details{These functions compute permutation adjusted \eqn{p}-values for the step-down maxT and minP multiple testing procedures, which provide strong control of the family-wise Type I error rate (FWER). The adjusted \eqn{p}-values for the minP procedure are defined in equation (2.10) p. 66 of Westfall & Young (1993), and the maxT procedure is discussed p. 50 and 114. The permutation algorithms for estimating the adjusted \eqn{p}-values are given in Ge et al. (In preparation). The procedures are for the simultaneous test of \eqn{m} null hypotheses, namely, the null hypotheses of no association between the \eqn{m} variables corresponding to the rows of the data frame \code{X} and the class labels \code{classlabel}. For gene expression data, the null hypotheses correspond to no differential gene expression across mRNA samples.
}


\value{
  A data frame with components
  \item{index}{Vector of row indices, between 1 and \code{nrow(X)}, where rows are sorted first according to
    their adjusted \eqn{p}-values, next their unadjusted \eqn{p}-values, and finally their test statistics. }
  \item{teststat}{Vector of test statistics, ordered according to \code{index}. To get the test statistics in the original data order, use \code{teststat[order(index)]}.}
  \item{rawp}{Vector of raw (unadjusted) \eqn{p}-values, ordered according to \code{index}.}
  \item{adjp}{Vector of adjusted \eqn{p}-values, ordered according to \code{index}.}
  \item{plower}{For \code{\link{mt.minP}} function only, vector of "adjusted \eqn{p}-values", where ties in the permutation distribution of the successive minima of raw \eqn{p}-values with the observed \eqn{p}-values are counted only once. Note that procedures based on \code{plower} do not control the FWER. Comparison of \code{plower} and \code{adjp} gives an idea of the discreteness of the permutation distribution. Values in \code{plower} are ordered according to \code{index}.}
}

\references{
S. Dudoit, J. P. Shaffer, and J. C. Boldrick (Submitted). Multiple hypothesis testing in microarray experiments.\cr

Y. Ge, S. Dudoit, and T. P. Speed. Resampling-based multiple testing for microarray data hypothesis, Technical Report \#633 of UCB Stat. \url{http://www.stat.berkeley.edu/~gyc} \cr

P. H. Westfall and S. S. Young (1993). \emph{Resampling-based
multiple testing: Examples and methods for \eqn{p}-value adjustment}. John Wiley \& Sons.
}	
  
\author{Yongchao Ge, \email{yongchao.ge@mssm.edu}, \cr
Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{mt.plot}}, \code{\link{mt.rawp2adjp}}, \code{\link{mt.reject}}, \code{\link{mt.sample.teststat}}, \code{\link{mt.teststat}}, \code{\link{golub}}.}


\examples{
# Gene expression data from Golub et al. (1999)
# To reduce computation time and for illustrative purposes, we condider only
# the first 100 genes and use the default of B=10,000 permutations.
# In general, one would need a much larger number of permutations
# for microarray data.

data(golub)
smallgd<-golub[1:100,] 
classlabel<-golub.cl

# Permutation unadjusted p-values and adjusted p-values 
# for maxT and minP procedures with Welch t-statistics
resT<-mt.maxT(smallgd,classlabel)
resP<-mt.minP(smallgd,classlabel)
rawp<-resT$rawp[order(resT$index)]
teststat<-resT$teststat[order(resT$index)]

# Plot results and compare to Bonferroni procedure
bonf<-mt.rawp2adjp(rawp, proc=c("Bonferroni"))
allp<-cbind(rawp, bonf$adjp[order(bonf$index),2], resT$adjp[order(resT$index)],resP$adjp[order(resP$index)])

mt.plot(allp, teststat, plottype="rvsa", proc=c("rawp","Bonferroni","maxT","minP"),leg=c(0.7,50),lty=1,col=1:4,lwd=2)
mt.plot(allp, teststat, plottype="pvsr", proc=c("rawp","Bonferroni","maxT","minP"),leg=c(60,0.2),lty=1,col=1:4,lwd=2)
mt.plot(allp, teststat, plottype="pvst", proc=c("rawp","Bonferroni","maxT","minP"),leg=c(-6,0.6),pch=16,col=1:4)

# Permutation adjusted p-values for minP procedure with F-statistics (like equal variance t-statistics)
mt.minP(smallgd,classlabel,test="f",fixed.seed.sampling="n")

# Note that the test statistics used in the examples below are not appropriate 
# for the Golub et al. data. The sole purpose of these examples is to 
# demonstrate the use of the mt.maxT and mt.minP functions.

# Permutation adjusted p-values for maxT procedure with paired t-statistics
classlabel<-rep(c(0,1),19)
mt.maxT(smallgd,classlabel,test="pairt")

# Permutation adjusted p-values for maxT procedure with block F-statistics
classlabel<-rep(0:18,2)
mt.maxT(smallgd,classlabel,test="blockf",side="upper")

}
\keyword{htest}
