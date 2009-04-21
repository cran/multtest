\name{EBMTP}
\alias{EBMTP}

\title{A function to perform empirical Bayes resampling-based multiple hypothesis testing}

\description{
A user-level function to perform empirical Bayes multiple testing procedures (EBMTP). A variety of t- and F-tests, including robust versions of most tests, are implemented.  A common-cutoff method is used to control the chosen type I error rate (FWER, gFWER, TPPFP, or FDR).  Bootstrap-based null distributions are available.  Additionally, for t-statistics, one may wish to sample from an appropriate multivariate normal distribution with mean zero and correlation matrix derived from the vector influence function.  In EBMTP, realizations of local q-values, obtained via density estimation, are used to partition null and observed test statistics into guessed sets of true and false null hypotheses at each round of (re)sampling.  In this manner, parameters of any type I error rate which can be expressed as a function the number of false positives and true positives can be estimated.  Arguments are provided for user control of output. Gene selection in microarray experiments is one application.
}

\usage{
EBMTP(X, W = NULL, Y = NULL, Z = NULL, Z.incl = NULL, Z.test = NULL, 
    na.rm = TRUE, test = "t.twosamp.unequalvar", robust = FALSE, 
    standardize = TRUE, alternative = "two.sided", typeone = "fwer", 
    method = "common.cutoff", k = 0, q = 0.1, alpha = 0.05, smooth.null = FALSE, 
    nulldist = "boot.cs", B = 1000, psi0 = 0, marg.null = NULL, 
    marg.par = NULL, ncp = NULL, perm.mat = NULL, ic.quant.trans = FALSE, 
    MVN.method = "mvrnorm", penalty = 1e-06, prior = "conservative", 
    bw = "nrd", kernel = "gaussian", seed = NULL, cluster = 1, 
    type = NULL, dispatch = NULL, keep.nulldist = TRUE, keep.rawdist = FALSE, 
    keep.falsepos = FALSE, keep.truepos = FALSE, keep.errormat = FALSE,
    keep.Hsets=FALSE, keep.margpar = TRUE, keep.index = FALSE, keep.label = FALSE) 
}

\arguments{
For brevity, the presentation of arguments below will highlight those which differ significantly from arguments in the other main-level user function \code{MTP}.  See \code{\link{MTP}} for further details.

  \item{typeone}{Character string indicating which type I error rate to control, by default family-wise error rate ('fwer'). Other options include generalized family-wise error rate ('gfwer'), with parameter \code{k} giving the allowed number of false positives, and tail probability of the proportion of false positives ('tppfp'), with parameter \code{q} giving the allowed proportion of false positives. The false discovery rate ('fdr') can also be controlled.  In particular, for 'gfwer', 'tppfp' and 'fdr', multiple testing is not performed via augmentation of the results of a FWER-controlling MTP.  Rather, using guessed sets of true and false null hypotheses, these error rates are controlled in a more direct manner.}
  \item{method}{Character string indicating the EBMTP method.  Currently only 'common.cutoff' is implemented.  This method is most similar to 'ss.maxT' in \code{MTP}.}
  \item{nulldist}{Character string indicating which resampling method to use for estimating the joint test statistics null distribution, by default the non-parametric bootstrap with centering and scaling ('boot.cs').  The old default 'boot' will still compile and will correspond to 'boot.cs'.  Other null distribution options include 'boot.ctr', 'boot.qt', and 'ic', corresponding to the centered-only bootstrap distribution, quantile-transformed bootstrap distribution, and influence curve multivariate normal joint null distribution, respectively.  The permutation distribution is not available.}
  \item{prior}{Character string indicating which choice of prior probability to use for estimating local q-values (i.e., the posterior probabilities of a null hypothesis being true given the value of its corresponding test statistic).  Default is 'conservative', in which case the prior is set to its most conservative value of 1, meaning that all hypotheses are assumed to belong to the set of true null hypotheses.  Other options include 'ABH' for the adaptive Benjamini-Hochberg estimator of the number/proportion of true null hypotheses, and 'EBLQV' for the empirical Bayes local q-value value estimator of the number/proportion of true null hypotheses.  If 'EBLQV', the estimator of the prior probability is taken to be the sum of the estimated local q-values divided by the number of tests.  Relaxing the prior may result in more rejections, albeit at a cost of type I error control under certain conditions.  See details and references.}
  \item{bw}{A character string argument to \code{density} indicating the smoothing bandwidth to be used during kernel density estimation.  Default is 'nrd'.}
  \item{kernel}{A character string argument to \code{density} specifying the smoothing kernel to be used.  Default is 'gaussian'.} 
  \item{keep.falsepos}{A logical indicating whether or not to store the matrix of guessed false positives at each round of (re)sampling.  The matrix has rows equal to the number of cut-offs (observed test statistics) and columns equal to the \code{B} number of bootstrap samples or samples from the multivariate normal distribution (if \code{nulldist='ic'}).  Default is 'FALSE'.}
  \item{keep.truepos}{A logical indicating whether or not to store the matrix of guessed true positives at each round of (re)sampling.  The matrix has rows equal to the number of cut-offs (observed test statistics) and columns equal to the \code{B} number of bootstrap samples or samples from the multivariate normal distribution (if \code{nulldist='ic'}).  Default is 'FALSE'.}
  \item{keep.errormat}{A logical indicating whether or not to store the matrix of type I error rate values at each round of (re)sampling.  The matrix has rows equal to the number of cut-offs (observed test statistics) and columns equal to the \code{B} number of bootstrap samples or samples from the multivariate normal distribution (if \code{nulldist='ic'}).  Default is 'FALSE'.  In the case of FDR-control, for example, this matrix is \code{falsepos}/(\code{falsepos} + \code{truepos}).  The row means of this matrix are eventually used for assigning/ordering adjusted p-values to test statistics of each hypothesis.}
   \item{keep.Hsets}{A logical indicating whether or not to return the matrix of indicators which partition the hypotheses into guessed sets of true and false null hypotheses at each round of (re)sampling.  Default is 'FALSE'.}
   \item{X, W, Y, Z, Z.incl, Z.test, na.rm, test, robust, standardize, alternative, k, q, alpha, smooth.null, B, psi0, marg.null, marg.par, ncp, perm.mat, ic.quant.trans, MVN.method, penalty, seed, cluster, type, dispatch, keep.nulldist, keep.rawdist, keep.margpar, keep.index, keep.label}{These arguments are all similarly used by the \code{MTP} function, and their use has been defined elsewhere.  Please consult the \code{link{MTP}} help file or the references for further details.  Note that the \code{MTP}-function arguments \code{get.cr, get.cutoff, get.adjp} are now DEPRECATED in the \code{EBMTP} function.  Only adjusted p-values are calculated by \code{EBMTP}. These adjusted p-values are returned in the same order as the original hypotheses and raw p-values (typically corresponding to rows of \code{X}.)} 
}

\details{
The EBMTP begins with a marginal nonparametric mixture model for estimating local q-values.  By definition, q-values are 'the opposite' of traditional p-values.  That is, q-values represent the probability of null hypothesis being true given the value of its corresponding test statistic.  If the test statistics Tn have marginal distribution f = pi*f_0 + (1-pi)f_1, where pi is the prior probability of a true null hypothesis and f_0 and f_1 represent the marginal null and alternative densities, respectively, then the local q-value function is given by pi*f_0(Tn)/f(Tn). \cr

One can estimate both the null density f_0 and full density f by applying kernel density estimation over the matrix of null test statistics and the vector of observed test statistics, respectively.  Practically, this step in \code{EBMTP} also ensures that sidedness is correctly accounted for among the test statistics and their estimated null distribution.  The prior probability pi can be set to its most conservative value of 1 or estimated by some other means, e.g., using the adaptive Benjamini Hochberg ('ABH') estimator or by summing up the estimated local q-values themselves ('EBLQV')and dividing by the number of tests.  Bounding these estimated probabilities by one provides a vector of estimated local q-values with length equal to the number of hypotheses.  Bernoulli 0/1 realizations of the posterior probabilities indicate which hypotheses are guessed as belonging to the true set of null hypotheses given the value of their test statistics.  Once this partitioning has been achieved, one can count the numbers of guessed false positives and guessed true positives at each round of (re)sampling that are obtained when using the value of an observed test statistic as a cut-off. \cr

EBMTPs use function closures to represent type I error rates in terms of their defining features.  Restricting the choice of type I error rate to 'fwer', 'gfwer', 'tppfp', and 'fdr', means that these features include whether to control the number of false positives or the proportion of false positives among the number of rejetions made (i.e., the false discovery proportion), whether we are controlling a tail probability or expected value error rate, and, in the case of tail probability error rates, what bound we are placing on the random variable defining the type I error rate (e.g., k for 'gfwer' or 'q' for 'tppfp').  Averaging the type I error results over B (bootstrap or multivariate normal) samples provides an estimator of the evidence against the null hypothesis (adjusted p-values) with respect to the choice of type I error rate.  Finally, a monotonicity constraint is placed on the adjusted p-values before being returned as output. \cr

As detailed in the references, relaxing the prior may result in a more powerful multiple testing procedure, albeit sometimes at the cost of type I error control.  Additionally, when the proportion of true null hypotheses is close to one, type I error control may also become an issue, even when using the most conservative prior probability of one.  This feature is known to occur with some other procedures which rely on the marginal nonparametric mixture model for estimating (local) q-values.  The slot \code{EB.h0M} returned by objects of class \code{EBMTP} is the sum of the local q-values estimated via kernel density estimation (divided by the total number of tests).  If this value is close to one (>0.9-0.95), the user will probably not want relax the prior, as even the conservative EBMTP might be approaching a performance bound with respect to type I error control.  The user is advised to begin by using the most 'conservative' prior, assess the estimated proportion of true null hypotheses, and then decide if relaxing the prior might be desired.  Gains in power over other multiple testing procedures have been observed even when using the most conservative prior of one. \cr  

Situations of moderate-high to high levels of correlation may also affect the results of multiple testing methods which use the same mixture model for generating q-values.  Microarray analysis represents a scenario in which dependence structures are typically weak enough to mitigate this concern.  On the other hand, the analysis of densely sampled SNPs, for example, may present problems.  \cr 
 
}

\value{
An object of class \code{EBMTP}.  Again, for brevity, the values below represent slots which distinguish objects of class \code{EBMTP} from those of class \code{MTP}. \cr

\item{\code{falsepos}}{A matrix with rows equal to the number of hypotheses and columns the number of samples of null test statistics (\code{B}) indicating the number of guessed false positives when using the corresponding value of the observed test statistic as a cut-off.  Not returned unless \code{keep.falsepos=TRUE}.}
\item{\code{truepos}}{A matrix with rows equal to the number of hypotheses and columns the number of samples of null test statistics (\code{B}) indicating the number of guessed true positives when using the corresponding value of the observed test statistic as a cut-off.  Not returned unless \code{keep.truepos=TRUE}.}
\item{\code{errormat}}{The matrix obtained after applying to type I error rate function closure to the matrices in \code{falsepos}, and, if applicable, \code{truepos}.  Not returned unless \code{keep.errormat=TRUE}.}
\item{\code{EB.h0M}}{The sum of the local q-values obtained after density estimation.  This number serves as an estimate of the proportion of true null hypotheses.  Values close to one indicate situations in which type I error control may not be guaranteed by the EBMTP.  When \code{prior='EBLQV'}, this value is used as the prior 'pi' during evaluation of the local q-value function.} 
\item{\code{prior}}{The numeric value of the prior 'pi' used when evaluating the local q-value function.}
\item{\code{prior.type}}{Character string returning the value of \code{prior} in the original call to \code{EBMTP}.  One of 'conservative', 'ABH', or 'EBLQV'.}
\item{\code{lqv}}{A numeric vector of length the number of hypotheses with the estimated local q-values used for generating guessed sets of true null hypotheses.}
\item{\code{Hsets}}{A numeric matrix with the same dimension as \code{nulldist}, containing the Bernoulli realizations of the estimated local q-values stored in \code{lqv} which were used to partition the hypotheses into guessed sets of true and false null hypotheses at each round of (re)sampling. Not returned unless \code{keep.Hsets=TRUE}.}
}

\references{

H.N. Gilbert, K.S. Pollard, M.J. van der Laan, and S. Dudoit (2009). Resampling-based multiple 
hypothesis testing with applications to genomics: New developments in R/Bioconductor 
package multtest. \emph{Journal of Statistical Software} (submitted). Temporary URL: \url{http://www.stat.berkeley.edu/~houston/JSSNullDistEBMTP.pdf}.\cr

Y. Benjamini and Y. Hochberg (2000). On the adaptive control of the false 
discovery rate in multiple testing with independent statistics. \emph{J. Behav.
Educ. Statist}. Vol 25: 60-83.\cr

Y. Benjamini, A. M. Krieger and D. Yekutieli (2006). Adaptive linear step-up
procedures that control the false discovery rate. \emph{Biometrika}. 
Vol. 93: 491-507.\cr

M.J. van der Laan, M.D. Birkner, and A.E. Hubbard (2005).  Empirical Bayes and Resampling Based Multiple Testing Procedure Controlling the Tail Probability of the Proportion of False Positives. Statistical Applications in Genetics and Molecular Biology, 4(1).
\url{http://www.bepress.com/sagmb/vol4/iss1/art29/} \cr

S. Dudoit and M.J. van der Laan.  Multiple Testing Procedures and Applications to Genomics.  Springer Series in Statistics. Springer, New York, 2008. \cr

S. Dudoit, H.N. Gilbert, and M J. van der Laan (2008). 
Resampling-based empirical Bayes multiple testing procedures for controlling 
generalized tail probability and expected value error rates: Focus on the false
discovery rate and simulation study. \emph{Biometrical Journal}, 50(5):716-44. \url{http://www.stat.berkeley.edu/~houston/BJMCPSupp/BJMCPSupp.html}. \cr

H.N. Gilbert, M.J. van der Laan, and S. Dudoit. Joint multiple testing procedures for 
graphical model selection with applications to biological networks. Technical report, 
U.C. Berkeley Division of Biostatistics Working Paper Series, April 2009. URL \url{http://www.bepress.com/ucbbiostat/paper245}. \cr


}

\author{Houston N. Gilbert, based on the original \code{MTP} code written by Katherine S. Pollard}

\seealso{\code{\link{MTP}}, \code{\link{EBMTP-class}}, \code{\link{EBMTP-methods}}, \code{\link{Hsets}}}

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
}
