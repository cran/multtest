\name{Hsets}
\alias{Hsets}
\alias{G.VS}
\alias{ABH.h0}
\alias{VScount}
\alias{dens.est}

\title{Functions for generating guessed sets of true null hypotheses in empirical Bayes resampling-based multiple hypothesis testing}

\description{These functions are called internally by the main user-level function \code{EBMTP}.  They are used for estimating local q-values, generating guessed sets of true null hypotheses, and applying these results to function closures defining the choice of type I error rate (FWER, gFWER, TPPFP, and FDR).}

\usage{
Hsets(Tn, nullmat, bw, kernel, prior, B, rawp) 

ABH.h0(rawp) 

G.VS(V, S = NULL, tp = TRUE, bound)
}

\arguments{
  \item{Tn}{The vector of observed test statistics.}
  \item{nullmat}{The matrix of null test statistics obtained either through null transformation of the bootstrap distribution or by sampling from an appropriate multivariate normal distribution (when \code{nulldist='ic'}.)}
  \item{bw}{A character string argument to \code{density} indicating the smoothing bandwidth to be used during kernel density estimation. Default is 'nrd'.} 
  \item{kernel}{A character string argument to \code{density} specifying the smoothing kernel to be used.  Default is 'gaussian'.} 
  \item{prior}{Character string indicating which choice of prior probability to use for estimating local q-values (i.e., the posterior probabilities of a null hypothesis being true given the value of its corresponding test statistic).  Default is 'conservative', in which case the prior is set to its most conservative value of 1, meaning that all hypotheses are assumed to belong to the set of true null hypotheses.  Other options include 'ABH' for the adaptive Benjamini-Hochberg estimator of the number/proportion of true null hypotheses, and 'EBLQV' for the empirical Bayes local q-value value estimator of the number/proportion of true null hypotheses.  If 'EBLQV', the estimator of the prior probability is taken to be the sum of the estimated local q-values divided by the number of tests.  Relaxing the prior may result in more rejections, albeit at a cost of type I error control under certain conditions.  See references.}
  \item{B}{The number of bootstrap iterations (i.e. how many resampled data sets) or the number of samples from the multivariate normal distribution (if \code{nulldist='ic'}). Can be reduced to increase the speed of computation, at a cost to precision. Default is 1000.}
  \item{rawp}{A vector of raw (unadjusted) p-values obtained bootstrap-based or influence curve null distribution.}
  \item{V}{A matrix of the numbers of guessed false positives for each cut-off, i.e., observed value of a test statistic, within each sample in \code{B}.}
  \item{S}{A matrix of the numbers of guessed true positives for each cut-off, i.e., observed value of a test statistic, within each sample in \code{B}.}
  \item{tp}{Logical indicator which is TRUE if type I error rate is a tail probability error rate and FALSE is if it is an expected value error rate.}
  \item{bound}{If a tail probability error rate, the bound to be placed on function of guessed false positives and guessed true positives.  For, 'fwer', equal to 0; 'gfwer', equal to 'k'; and tppfp, equal to 'q'.} 
}

\details{
The most important object to be returned from the function \code{Hsets} is a matrix of indicators, i.e., Bernoulli realizations of the estimated local q-values, taking the value of 1 if the hypothesis is guessed as belonging to the set of true null hypotheses and 0 otherwise (guessed true alternative).  Realizations of these probabilities are generated with a call to \code{rbinom}, meaning that this function will set the RNG seed forward another \code{B}*(the number of hypotheses) places.  This matrix, with rows equal to the number of hypotheses and columns the number of (bootstrap or multivariate normal) samples is used to subset the matrix of null test statistics and the vector of observed test statistics at each round of (re)sampling into samples of statistics guessed as belonging to the sets of true null and true alternative hypotheses, respectively.  Using the values of the observed test statistics themselves as cut-offs, the numbers of guessed false positives and (if applicable) guessed true positives can be counted and eventually used to estimate the distribution of a type I error rate characterized by the closure returned from \code{G.VS}.  Counting of guessed false positives and guessed true positives is performed in C through the function \code{VScount}.  
}

\value{
For the function \code{Hsets}, a list with the following elements:
  \item{Hsets.mat}{A matrix of numeric indicators with rows equal to the number of test (hypotheses, typically \code{nrow(X)}) and columns the number of samples of null test statistics, \code{B}.  Values of one indicate hypotheses guessed as belonging to the set of true null hypotheses based on the value of their corresponding test statistic.  Values of zero correspond to hypotheses guesses as belonging to the set of true alternative hypotheses.}
  \item{EB.h0M}{The estimated proportion of true null hypotheses as determined by nonparametric density estimation.  This value is the sum of the estimated local q-values divided by the total number of tests (hypotheses).} 
  \item{prior}{The value of the prior applied to the local q-value function.  If 'conservative', the prior is set to one.  Otherwise, the prior is the value obtained from the estimator of the adaptive Benjamini-Hochberg procedure (if \code{prior} is 'ABH') or from density estimation (if \code{prior} is 'EBLQV').}
  \item{pn.out}{The vector of estimated local q-values.  This vector is returned in the \code{lqv} slot of objects of class \code{EBMTP}.}

For the function \code{ABH.h0}, the estimated number of true null hypotheses using the estimator from the linear step-up adaptive Benjamini-Hochberg procedure. \cr

For the function \code{G.VS}, a closure which accepts as arguments the matrices of guessed false positive and true positives (if applicable) and applies the appropriate function defining the desired type I error rate.
}

\references{
H.N. Gilbert, K.S. Pollard, M.J. van der Laan, and S. Dudoit (2009). Resampling-based multiple 
hypothesis testing with applications to genomics: New developments in R/Bioconductor 
package multtest. \emph{Journal of Statistical Software} (submitted). Temporary URL: \url{http://www.stat.berkeley.edu/~houston/JSSNullDistEBMTP.pdf}.\cr

Y. Benjamini and Y. Hochberg (2000). On the adaptive control of the false 
discovery rate in multiple testing with independent statistics. \emph{J. Behav.
Educ. Statist}. Vol 25: 60-83.\cr

Y. Benjamini, A.M. Krieger and D. Yekutieli (2006). Adaptive linear step-up
procedures that control the false discovery rate. \emph{Biometrika}. 
Vol. 93: 491-507.\cr

M.J. van der Laan, M.D. Birkner, and A.E. Hubbard (2005).  Empirical Bayes and Resampling Based Multiple Testing Procedure Controlling the Tail Probability of the Proportion of False Positives. Statistical Applications in Genetics and Molecular Biology, 4(1).
\url{http://www.bepress.com/sagmb/vol4/iss1/art29/} \cr

S. Dudoit and M.J. van der Laan.  Multiple Testing Procedures and Applications to Genomics.  Springer Series in Statistics. Springer, New York, 2008. \cr

S. Dudoit, H.N. Gilbert, and M.J. van der Laan (2008). 
Resampling-based empirical Bayes multiple testing procedures for controlling 
generalized tail probability and expected value error rates: Focus on the false
discovery rate and simulation study. \emph{Biometrical Journal}, 50(5):716-44. \url{http://www.stat.berkeley.edu/~houston/BJMCPSupp/BJMCPSupp.html}. \cr

H.N. Gilbert, M.J. van der Laan, and S. Dudoit. Joint multiple testing procedures for 
graphical model selection with applications to biological networks. Technical report, 
U.C. Berkeley Division of Biostatistics Working Paper Series, April 2009. URL \url{http://www.bepress.com/ucbbiostat/paper245}. \cr
}

\author{Houston N. Gilbert}

\seealso{\code{\link{EBMTP}}, \code{\link{EBMTP-class}}, \code{\link{EBMTP-methods}}}

\examples{
set.seed(99)
data<-matrix(rnorm(90),nr=9)
group<-c(rep(1,5),rep(0,5))

#EB fwer control with centered and scaled bootstrap null distribution 
#(B=100 for speed)
eb.m1<-EBMTP(X=data,Y=group,alternative="less",B=100,method="common.cutoff")
print(eb.m1)
summary(eb.m1)
par(mfrow=c(2,2))
plot(eb.m1,top=9)

abh <- ABH.h0(eb.m1@rawp)
abh

eb.m2 <- EBupdate(eb.m1,prior="ABH")
eb.m2@prior

}


