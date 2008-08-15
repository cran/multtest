\name{MTP}
\alias{MTP}

\title{A function to perform resampling-based multiple hypothesis testing}

\description{
A user-level function to perform multiple testing procedures (MTP). A variety of t- and f-tests, including robust versions of each test, are implemented. Single-step and step-down minP and maxT methods are used to control the chosen type I error rate (FWER, gFWER, TPPFP, or FDR). Bootstrap and permutation null distributions are available. Arguments are provided for user control of output. Gene selection in microarray experiments is one application.
}

\usage{
MTP(X, W = NULL, Y = NULL, Z = NULL, Z.incl = NULL, Z.test = NULL, 
na.rm = TRUE, test = "t.twosamp.unequalvar", robust = FALSE, 
standardize = TRUE, alternative = "two.sided", psi0 = 0, typeone = "fwer", 
k = 0, q = 0.1, fdr.method = "conservative", alpha = 0.05, smooth.null = 
FALSE, nulldist = "boot", csnull=TRUE,B = 1000, method = "ss.maxT", get.cr = FALSE, 
get.cutoff = FALSE, get.adjp = TRUE, keep.nulldist = TRUE, seed = NULL, cluster=1,
type = NULL, dispatch = NULL)
}

\arguments{
  \item{X}{A matrix, data.frame or ExpressionSet containing the raw data. In the case of an ExpressionSet, \code{exprs(X)} is the data of interest and \code{pData(X)} may contain outcomes and covariates of interest. For currently implemented tests, one hypothesis is tested for each row of the data.}
  \item{W}{A vector or matrix containing non-negative weights to be used in computing the test statistics. If a matrix, \code{W} must be the same dimension as \code{X} with one weight for each value in \code{X}. If a vector, \code{W} may contain one weight for each observation (i.e. column) of \code{X} or one weight for each variable (i.e. row) of \code{X}. In either case, the weights are duplicated appropraiately. Weighted f-tests are not available. Default is 'NULL'.}
  \item{Y}{A vector, factor, or \code{Surv} object containing the outcome of interest. This may be class labels (f-tests and two sample t-tests) or a continuous or polycotomous dependent variable (linear regression based t-tests), or survival data (Cox proportional hazards based t-tests). For \code{block.f} and \code{f.twoway} tests, class labels must be ordered by block and within each block ordered by group. If \code{X} is an ExpressionSet, \code{Y} can be a character string referring to the column of \code{pData(X)} to use as outcome. Default is 'NULL'.}
  \item{Z}{A vector, factor, or matrix containing covariate data to be used in the regression (linear and Cox) models. Each variable should be in one column, so that \code{nrow(Z)=ncol(X)}. If \code{X} is an ExpressionSet, \code{Z} can be a character string referring to the column of \code{pData(X)} to use as covariates. The variables \code{Z.incl} and \code{Z.adj} allow one to specify which covariates to use in a particular test without modifying the input \code{Z}. Default is 'NULL'.}
  \item{Z.incl}{The indices of the columns of \code{Z} (i.e. which variables) to include in the model. These can be numbers or column names (if the columns are names). Default is 'NULL'.}
  \item{Z.test}{The index or names of the column of \code{Z} (i.e. which variable) to use to test for association with each row of \code{X} in a linear model. Only used for \code{test="lm.XvsZ"}, where it is necessary to specify which covariate's regression parameter is of interest. Default is 'NULL'.}
  \item{na.rm}{Logical indicating whether to remove observations with an NA. Default is 'TRUE'.}
  \item{test}{Character string specifying the test statistics to use, by default 't.twosamp.unequalvar'. See details (below) for a list of tests.}
  \item{robust}{Logical indicating whether to use the robust version of the chosen test, e.g. Wilcoxon singed rank test for robust one-sample t-test or \code{rlm} instead of \code{lm} in linear models. Default is 'FALSE'.}
  \item{standardize}{Logical indicating whether to use the standardized version of the test statistics (usual t-statistics are standardized). Default is 'TRUE'.}
  \item{alternative}{Character string indicating the alternative hypotheses, by default 'two.sided'. For one-sided tests, use 'less' or 'greater' for null hypotheses of 'greater than or equal' (i.e. alternative is 'less') and 'less than or equal', respectively.}
  \item{psi0}{The hypothesized null value, typically zero (default). Currently, this should be a single value, which is used for all hypotheses.}
  \item{typeone}{Character string indicating which type I error rate to control, by default family-wise error rate ('fwer'). Other options include generalized family-wise error rate ('gfwer'), with parameter \code{k} giving the allowed number of false positives, and tail probability of the proportion of false positives ('tppfp'), with parameter \code{q} giving the allowed proportion of false positives. The false discovery rate ('fdr') can also be conrtolled.}
  \item{k}{The allowed number of false positives for gFWER control. Default is 0 (FWER).}
  \item{q}{The allowed proportion of false positives for TPPFP control. Default is 0.1.}
  \item{fdr.method}{Character string indicating which FDR controlling method should be used when \code{typeone="fdr"}. The options are "conservative" (default) for the more conservative, general FDR controlling procedure and "restricted" for the method which requires more assumptions.}
  \item{alpha}{The target nominal type I error rate, which may be a vector of error rates. Default is 0.05.}
  \item{smooth.null}{Indicator of whether to use a kernal density estimate for the tail of the null distributon for computing raw pvalues close to zero. Only used if 'rawp' would be zero without smoothing. Default is 'FALSE'.}
  \item{nulldist}{Character string indicating which resampling method to use for estimating the joint test statistics null distribution, by default non-parametric bootstrap ('boot').}
  \item{csnull}{Indicator of whether the bootstrap estimated test statistics distribution should be centered and scaled (to produce a null distirbution) or not. If \code{csnull==FALSE}, the non-null bootstrap estimated test statistics distribution is returned.} 
  \item{B}{The number of bootstrap iterations (i.e. how many resampled data sets) or the number of permutations (if \code{nulldist} is 'perm'). Can be reduced to increase the speed of computation, at a cost to precision. Default is 1000.}
  \item{method}{The multiple testing procedure to use. Options are single-step maxT ('ss.maxT', default), single-step minP ('ss.minP'), step-down maxT ('sd.maxT'), and step-down minP ('sd.minP').}
  \item{get.cr}{Logical indicating whether to compute confidence intervals for the estimates. Not available for f-tests. Default is 'FALSE'.}
  \item{get.cutoff}{Logical indicating whether to compute thresholds for the test statistics. Default is 'FALSE'.}
  \item{get.adjp}{Logical indicating whether to compute adjusted p-values. Default is 'TRUE'.}
  \item{keep.nulldist}{Logical indicating whether to return the computed null distribution, by default 'TRUE'. Note that this matrix can be quite large. }
  \item{seed}{Integer or vector of integers to be used as argument to \code{set.seed} to set the seed for the random number generator for bootstrap resampling. This argument can be used to repeat exactly a test performed with a given seed. If the seed is specified via this argument, the same seed will be returned in the seed slot of the MTP object created. Else a random seed(s) will be generated, used and returned. Vector of integers used to specify seeds for each node in a cluster used to to generate a bootstrap null distribution.}
  \item{cluster}{Integer for number of nodes to create or a cluster object created through the package snow. With \code{cluster=1}, bootstrap is implemented on single node. Supplying a cluster object results in the bootstrap being implemented in parallel on the provided nodes. This option is only available for the bootstrap procedure. With default value of 1, bootstrap is executed on single CPU.}
  \item{type}{Interface system to use for computer cluster. See \code{snow} package for details.}
  \item{dispatch}{The number or percentage of bootstrap iterations to dispatch at a time to each node of the cluster if a computer cluster is used. If dispatch is a percentage, \code{B*dispatch} must be an integer. If dispatch is an integer, then \code{B/dispatch} must be an integer. Default is 5 percent.}
}

\details{
A multiple testing procedure (MTP) is defined by choices of test statistics, type I error rate, null distribution and method for error rate control. Each component is described here. See references for more detail.

Test statistics are determined by the values of \code{test}: 
\describe{
\item{t.onesamp:}{one-sample t-statistic for tests of means;}
\item{t.twosamp.equalvar:}{equal variance two-sample t-statistic for tests of differences in means (two-sample t-statistic);}
\item{t.twosamp.unequalvar:}{unequal variance two-sample t-statistic for tests of differences in means (two-sample Welch t-statistic);}
\item{t.pair:}{two-sample paired t-statistic for tests of differences in means;}
\item{f:}{multi-sample f-statistic for tests of equality of population means (assumes constant variance across groups, but not normality); }
\item{f.block:}{multi-sample f-statistic for tests of equality of population means in a block design (assumes constant variance across groups, but not normality). This test is not available with the bootstrap null distribution;}
\item{f.twoway:}{multi-sample f-statistic for tests of equality of population means in a block design (assumes constant variance across groups, but not normality). Differs from \code{f.block} in requiring multiple observations per group*block combintation. This test uses the means of each group*block combination as response variable and test for group main effects assuming a randomized block design;}
\item{lm.XvsZ:}{t-statistic for tests of regression coefficients for variable \code{Z.test} in linear models, each with a row of X as outcome, possibly adjusted by covariates \code{Z.incl} from the matrix \code{Z} (in the case of no covariates, one recovers the one-sample t-statistic, \code{t.onesamp});}
\item{lm.YvsXZ:}{t-statistic for tests of regression coefficients in linear models, with outcome Y and each row of X as covariate of interest, with possibly other covariates \code{Z.incl} from the matrix \code{Z};}
\item{coxph.YvsXZ:}{t-statistic for tests of regression coefficients in Cox proportional hazards survival models, with outcome Y and each row of X as covariate of interest, with possibly other covariates \code{Z.incl} from the matrix \code{Z}.}
}

When \code{robust=TRUE}, non-parametric versions of each test are performed. For the linear models, this means \code{rlm} is used instead of \code{lm}. There is not currently a robust version of \code{test=coxph.YvsXZ}. For the t- and f-tests, data values are simply replaced by their ranks. This is equivalent to performing the following familiar named rank-based tests. The conversion after each test is the formula to convert from the MTP test to the statistic reported by the listed R function (where num is the numerator of the MTP test statistics, n is total sample size, nk is group k sample size, K is total number of groups or treatments, and rk are the ranks in group k).
\describe{
\item{t.onesamp or t.pair:}{Wilcoxon signed rank, \code{wilcox.test} with \code{y=NULL} or \code{paired=TRUE}, \cr
conversion: num/n}
\item{t.twosamp.equalvar:}{Wilcoxon rank sum or Mann-Whitney, \code{wilcox.test}, \cr
conversion: n2*(num+mean(r1)) - n2*(n2+1)/2}
\item{f:}{Kruskal-Wallis rank sum, \code{kruskal.test}, \cr
conversion: num*12/(n*(n-1)}
\item{f.block:}{Friedman rank sum, \code{friedman.test}, \cr
conversion: num*12/(K*(K+1))}
\item{f.twoway:}{Friedman rank sum, \code{friedman.test}, \cr
conversion: num*12/(K*(K+1))}
}

The implemented MTPs are based on control of the family-wise error rate, defined as the probability of any false positives. Let Vn denote the (unobserved) number of false positives. Then, control of FWER at level alpha means that Pr(Vn>0)<=alpha. The set of rejected hypotheses under a FWER controlling procedure can be augmented to increase the number of rejections, while controlling other error rates. The generalized family-wise error rate is defined as Pr(Vn>k)<=alpha, and it is clear that one can simply take an FWER controlling procedure, reject k more hypotheses and have control of gFWER at level alpha. The tail probability of the proportion of false positives depends on both the number of false postives (Vn) and the number of rejections (Rn). Control of TPPFP at level alpha means Pr(Vn/Rn>q)<=alpha, for some proportion q. Control of the false discovery rate refers to the expected proportion of false positives (rather than a tail probability). Control of FDR at level alpha means E(Vn/Rn)<=alpha.

In practice, one must choose a method for estimating the test statistics null distribution. We have implemented an ordinary non-parametric bootstrap estimator and a permutation estimator (which makes sense in certain settings, see references). The non-parametric bootstrap estimator (default) provides asymptotic control of the type I error rate for any data generating distribution, whereas the permutation estimator requires the subset pivotality assumption. One draw back of both methods is the discreteness of the estimated null distribution when the sample size is small. Furthermore, when the sample size is small enough, it is possible that ties will lead to a very small variance estimate. Using \code{standardize=FALSE} allows one to avoid these unusually small test statistic denominators. Parametric bootstrap estimators are another option (not yet implemented).

Given observed test statistics, a type I error rate (with nominal level), and a test statistics null distribution, MTPs provide adjusted p-values, cutoffs for test statistics, and possibly confidence regions for estimates. Four methods are implemented, based on minima of p-values and maxima of test statistics. Only the step down methods are currently available with the permutation null distribution.

Computation times using a bootstrap null distribution are slower when weights are used for one and two-sample tests. Computation times when using a bootstrap null distribution also are slower for the tests \code{lmXvsZ}, \code{lmYvsXZ}, \code{coxph.YvsXZ}.

To execute the bootstrap on a computer cluster, a cluster object generated with \code{makeCluster} in the package \code{snow} may be used as the argument for cluster. Alternatively, the number of nodes to use in the computer cluster can be used as the argument to cluster. In this case, \code{type} must be specified and a cluster will be created. In both cases, \code{Biobase} and \code{multtest} will be loaded onto each cluster node if these libraries are located in a directory in the standard search path. If these libraries are in a non-standard location, it is necessary to first create the cluster, load \code{Biobase} and \code{multtest} on each node and then to use the cluster object as the argument to cluster. See documentation for \code{snow} package for additional information on creating and using a cluster.

}

\value{
An object of class \code{MTP}, with the following slots:

    \item{\code{statistic}}{Object of class \code{numeric}, observed test statistics for each hypothesis, specified by the values of the \code{MTP} arguments \code{test}, \code{robust}, \code{standardize}, and \code{psi0}.}
    \item{\code{estimate}}{For the test of single-parameter null hypotheses using t-statistics (i.e., not the F-tests), the numeric vector of estimated parameters corresponding to each hypothesis, e.g. means, differences in means, regression parameters.}
    \item{\code{sampsize}}{Object of class \code{numeric}, number of columns (i.e. observations) in the input data set.}
    \item{\code{rawp}}{Object of class \code{numeric}, unadjusted, marginal p-values for each hypothesis.}
    \item{\code{adjp}}{Object of class \code{numeric}, adjusted (for multiple testing) p-values for each hypothesis (computed only if the \code{get.adjp} argument is TRUE).}
    \item{\code{conf.reg}}{For the test of single-parameter null hypotheses using t-statistics (i.e., not the F-tests), the numeric array of lower and upper simultaneous confidence limits for the parameter vector, for each value of the nominal Type I error rate \code{alpha} (computed only if the \code{get.cr} argument is TRUE).}
    \item{\code{cutoff}}{The numeric matrix of cut-offs for the vector of test statistics for each value of the nominal Type I error rate \code{alpha} (computed only if the \code{get.cutoff} argument is TRUE).}
   \item{\code{reject}}{Object of class \code{"matrix"}, rejection indicators (TRUE for a rejected null hypothesis), for each value of the nominal Type I error rate \code{alpha}.}
    \item{\code{nulldist}}{The numeric matrix for the estimated test statistics null distribution (returned only if \code{keep.nulldist=TRUE}; option not currently available for permutation null distribution, i.e.,  \code{nulldist="perm"}). By default (i.e., for \code{nulldist="boot"}), the entries of \code{nulldist} are the null value shifted and scaled bootstrap test statistics, with one null test statistic value for each hypothesis (rows) and bootstrap iteration (columns).}
   \item{\code{call}}{Object of class \code{call}, the call to the MTP function.}
    \item{\code{seed}}{An integer or vector for specifying the state of the random number generator used to create the resampled datasets. The seed can be reused for reproducibility in a repeat call to \code{MTP}. This argument is currently used only for the bootstrap null distribution (i.e., for \code{nulldist="boot"}). See \code{? set.seed} for details.}
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

\author{Katherine S. Pollard with design contributions from Sandra Taylor, Sandrine Dudoit and Mark J. van der Laan.}

\note{Thank you to Peter Dimitrov for suggestions about the code.}

\seealso{\code{\link{MTP-class}}, \code{\link{MTP-methods}}, \code{\link{mt.minP}}, \code{\link{mt.maxT}}, \code{\link{ss.maxT}}, \code{\link{fwer2gfwer}}}

\examples{

#data 
set.seed(99)
data<-matrix(rnorm(90),nr=9)
group<-c(rep(1,5),rep(0,5))

#fwer control with bootstrap null distribution (B=100 for speed)
m1<-MTP(X=data,Y=group,alternative="less",B=100,method="sd.minP")
print(m1)
summary(m1)
par(mfrow=c(2,2))
plot(m1,top=9)

}

\keyword{htest}


