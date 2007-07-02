\name{mt.plot}
\alias{mt.plot}
\title{Plotting results from multiple testing procedures}
\description{This function produces a number of graphical summaries
  for the results of multiple testing procedures and their corresponding
  adjusted \eqn{p}-values.}
\usage{
mt.plot(adjp, teststat, plottype="rvsa", logscale=FALSE, alpha=seq(0, 1, length = 100), proc, leg=c(0, 0), \dots)
}
\arguments{
  \item{adjp}{A matrix of adjusted \emph{p}-values, with rows
    corresponding to hypotheses (genes) and columns to multiple testing
    procedures. This matrix could be obtained from the functions
    \code{\link{mt.maxT}}, \code{\link{mt.minP}}, or \code{\link{mt.rawp2adjp}}.}
  \item{teststat}{A vector of test statistics for each of the hypotheses. This vector could be obtained from the functions \code{\link{mt.teststat}}, \code{\link{mt.maxT}}, or \code{\link{mt.minP}}.}
  \item{plottype}{A character string specifying the type of graphical
    summary for the results of the multiple testing procedures. \cr
    If \code{plottype="rvsa"}, the number of rejected hypotheses is plotted against the nominal Type I error rate for each of the procedures given in \code{proc}.\cr
    If \code{plottype="pvsr"}, the ordered adjusted \emph{p}-values are plotted for each of the procedures given in \code{proc}. This can be viewed as a plot of the Type I error rate against the number of rejected hypotheses. \cr
    If \code{plottype="pvst"}, the adjusted \emph{p}-values are plotted against the test statistics for each of the procedures given in \code{proc}.
    \cr
    If \code{plottype="pvsi"}, the adjusted \emph{p}-values are plotted for each of the procedures given in \code{proc} using the original data order.  }
  \item{logscale}{A logical variable for the \code{pvst} and \code{pvsi} plots. If \code{logscale} is \code{TRUE}, the negative decimal logarithms of the adjusted \emph{p}-values are plotted against the test statistics or gene indices. If \code{logscale} is \code{FALSE}, the adjusted \emph{p}-values are plotted against the test statistics or gene indices.}
  \item{alpha}{A vector of nominal Type I error rates for the \code{rvsa} plot.}
  \item{proc}{A vector of character strings containing the names of the
    multiple testing procedures, to be used in the legend.}
  \item{\dots}{Graphical parameters such as \code{col}, \code{lty},
    \code{pch}, and \code{lwd}
    may also be supplied as arguments to the function (see \code{\link{par}}).}
  \item{leg}{A vector of coordinates for the legend.}
}

\references{
  
  S. Dudoit, J. P. Shaffer, and J. C. Boldrick (Submitted). Multiple hypothesis testing in microarray experiments.\cr
  
  Y. Ge, S. Dudoit, and T. P. Speed. Resampling-based multiple testing for microarray data hypothesis, Technical Report \#633 of UCB Stat. \url{http://www.stat.berkeley.edu/~gyc}  \cr

}
\author{
  Sandrine Dudoit,  \url{http://www.stat.berkeley.edu/~sandrine}, \cr
  Yongchao Ge, \email{yongchao.ge@mssm.edu}.
}
\seealso{\code{\link{mt.maxT}}, \code{\link{mt.minP}}, \code{\link{mt.rawp2adjp}},  \code{\link{mt.reject}}, \code{\link{mt.teststat}}, \code{\link{golub}}.}

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
res1<-mt.maxT(smallgd,classlabel)
rawp<-res1$rawp[order(res1$index)]
teststat<-res1$teststat[order(res1$index)]

# Permutation adjusted p-values for simple multiple testing procedures
procs<-c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY")
res2<-mt.rawp2adjp(rawp,procs)

# Plot results from all multiple testing procedures
allp<-cbind(res2$adjp[order(res2$index),],res1$adjp[order(res1$index)])
dimnames(allp)[[2]][9]<-"maxT"
procs<-dimnames(allp)[[2]]
procs[7:9]<-c("maxT","BH","BY")
allp<-allp[,procs]

cols<-c(1:4,"orange","brown","purple",5:6)
ltypes<-c(3,rep(1,6),rep(2,2))

# Ordered adjusted p-values
mt.plot(allp,teststat,plottype="pvsr",proc=procs,leg=c(80,0.4),lty=ltypes,col=cols,lwd=2)

# Adjusted p-values in original data order
mt.plot(allp,teststat,plottype="pvsi",proc=procs,leg=c(80,0.4),lty=ltypes,col=cols,lwd=2)

# Number of rejected hypotheses vs. level of the test
mt.plot(allp,teststat,plottype="rvsa",proc=procs,leg=c(0.05,100),lty=ltypes,col=cols,lwd=2)

# Adjusted p-values vs. test statistics
mt.plot(allp,teststat,plottype="pvst",logscale=TRUE,proc=procs,leg=c(0,4),pch=ltypes,col=cols)

}
\keyword{hplot}
