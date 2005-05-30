\name{ss.maxT}
\alias{ss.maxT}
\alias{ss.minP}
\alias{sd.maxT}
\alias{sd.minP}

\title{Procedures to perform multiple testing}

\description{
Given observed test statistics, a test statistics null distribution, and alternetive hyptheses, these multiple testing procedures provide family-wise error rate (FWER) adjusted p-values, cutoffs for test statistics, and possibly confidence regions for estimates. Four methods are implemented, based on minima of p-values and maxima of test statistics. 
}

\usage{
ss.maxT(null, obs, alternative, get.cutoff, get.cr, 
get.adjp, alpha = 0.05)

ss.minP(null, obs, rawp, alternative, get.cutoff, get.cr, 
get.adjp, alpha=0.05)

sd.maxT(null, obs, alternative, get.cutoff, get.cr, 
get.adjp, alpha = 0.05)

sd.minP(null, obs, rawp, alternative, get.cutoff, get.cr, 
get.adjp, alpha=0.05)
}

\arguments{
  \item{null}{A matrix containing the test statistics null distribution, e.g. the output of \code{boot.resample}.}
  \item{obs}{A vector of observed test statistics, e.g. the output of a test statistics closure such as \code{meanX}. These are stored as a matrix with numerator (possibly absolute value or negative, depending on the value of alternative) in the first row, denominator in the second row, and a 1 or -1 in the third row (depending on the value of alternative). The observed test statistics are obs[1,]*obs[3,]/obs[2,].}
  \item{rawp}{Numeric vector of unadjusted ("raw") marginal p-values.}
  \item{alternative}{Character string indicating the alternative hypotheses, by default 'two.sided'. For one-sided tests, use 'less' or 'greater' for null hypotheses of 'greater than or equal' (i.e. alternative is 'less') and 'less than or equal', respectively.}
  \item{get.cutoff}{Logical indicating whether to compute thresholds for the test statistics. Default is 'FALSE'.}
  \item{get.cr}{Logical indicating whether to compute confidence intervals for the estimates. Not available for f-tests. Default is 'FALSE'.}
  \item{get.adjp}{Logical indicating whether to compute adjusted p-values. Default is 'TRUE'.}
  \item{alpha}{The target nominal type I error rate, which may be a vector of error rates. Default is 0.05.}
}

\details{
Having selected a suitable test statistics null distribution, there remains the main task of specifying rejection regions for each null hypothesis, i.e., cut-offs for each test statistic. One usually distinguishes between two main classes of multiple testing procedures, single-step and stepwise procedures. In single-step procedures, each null hypothesis is evaluated using a rejection region that is  independent of the results of the tests of other hypotheses. Improvement in power, while preserving Type I error rate control, may be achieved by stepwise (step-down or step-up) procedures, in which rejection of a particular null hypothesis depends on the outcome of the tests of other hypotheses. That is, the (single-step) test procedure is applied to a sequence of successively smaller nested random (i.e., data-dependent) subsets of null hypotheses, defined by the ordering of the test statistics (common cut-offs or maxT procedures) or unadjusted p-values (common-quantiles or minP procedures). 

In step-down procedures, the hypotheses corresponding to the most significant test statistics (i.e., largest absolute test statistics or smallest unadjusted p-values) are considered successively, with further tests depending on the outcome of earlier ones. As soon as one fails to reject a null hypothesis, no further
hypotheses are rejected. In contrast, for step-up procedures, the hypotheses corresponding to the least significant test statistics are considered successively, again with further tests depending on the outcome of earlier ones. As soon as one hypothesis is rejected, all remaining more significant hypotheses are rejected.

These functions perform the following procedures: \cr
ss.maxT: single-step, common cut-off (maxima of test statistics) \cr
ss.minP: single-step, common quantile (minima of p-values) \cr
sd.maxT: step-down, common cut-off (maxima of test statistics) \cr
sd.minP: step-down, common quantile (minima of p-values) \cr

}

\value{A list with the following components:
  \item{c}{Object of class \code{"matrix"}, for each nominal (i.e. target) level for the test, a vector of threshold values for the vector of test statistics.}
  \item{cr}{Object of class \code{"array"}, for each nominal (i.e. target) level for the test, a matrix of lower and upper confidence bounds for the parameter of interest for each hypothesis. Not available for f-tests.}
  \item{adjp}{Object of class \code{"numeric"}, adjusted p-values for each hypothesis.}
}

\references{
M.J. van der Laan, S. Dudoit, K.S. Pollard (2004), Augmentation Procedures for Control of the Generalized Family-Wise Error Rate and Tail Probabilities for the Proportion of False Positives, Statistical Applications in Genetics and Molecular Biology, 3(1). 
\url{http://www.bepress.com/sagmb/vol3/iss1/art15/}

M.J. van der Laan, S. Dudoit, K.S. Pollard (2004), Multiple Testing. Part II. Step-Down Procedures for Control of the Family-Wise Error Rate, Statistical Applications in Genetics and Molecular Biology, 3(1).
\url{http://www.bepress.com/sagmb/vol3/iss1/art14/}

S. Dudoit, M.J. van der Laan, K.S. Pollard (2004), Multiple Testing. Part I. Single-Step Procedures for Control of General Type I Error Rates, Statistical Applications in Genetics and Molecular Biology, 3(1).
\url{http://www.bepress.com/sagmb/vol3/iss1/art13/}

Katherine S. Pollard and Mark J. van der Laan, "Resampling-based Multiple Testing: Asymptotic Control of Type I Error and Applications to Gene Expression Data" (June 24, 2003). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 121.
\url{http://www.bepress.com/ucbbiostat/paper121}
}

\author{Katherine S. Pollard, \url{http://lowelab.ucsc.edu/katie/}\cr
with design contributions from Sandrine Dudoit and Mark J. van der Laan.}

\seealso{\code{\link{MTP}}}

\examples{
## These functions are used internally by the MTP function
## See MTP function: ? MTP
}

\keyword{htest}
\keyword{internal}



