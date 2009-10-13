\name{MTP}
\alias{MTP}

\title{A function to perform resampling-based multiple hypothesis testing}

\description{
A user-level function to perform multiple testing procedures (MTP). A variety of t- and F-tests, including robust versions of most tests, are implemented. Single-step and step-down minP and maxT methods are used to control the chosen type I error rate (FWER, gFWER, TPPFP, or FDR). Bootstrap and permutation null distributions are available.  Additionally, for t-statistics, one may wish to sample from an appropriate multivariate normal distribution with mean zero and correlation matrix derived from the vector influence function.  Arguments are provided for user control of output. Gene selection in microarray experiments is one application.
}

\usage{
MTP(X, W = NULL, Y = NULL, Z = NULL, Z.incl = NULL, Z.test = NULL, 
    na.rm = TRUE, test = "t.twosamp.unequalvar", robust = FALSE, 
    standardize = TRUE, alternative = "two.sided", psi0 = 0, 
    typeone = "fwer", k = 0, q = 0.1, fdr.method = "conservative", 
    alpha = 0.05, smooth.null = FALSE, nulldist = "boot.cs", 
    B = 1000, ic.quant.trans = FALSE, MVN.method = "mvrnorm", 
    penalty = 1e-06, method = "ss.maxT", get.cr = FALSE, get.cutoff = FALSE, 
    get.adjp = TRUE, keep.nulldist = TRUE, keep.rawdist = FALSE, 
    seed = NULL, cluster = 1, type = NULL, dispatch = NULL, marg.null = NULL, 
    marg.par = NULL, keep.margpar = TRUE, ncp = NULL, perm.mat = NULL, 
    keep.index = FALSE, keep.label = FALSE) 
}

\arguments{
  \item{X}{A matrix, data.frame or ExpressionSet containing the raw data. In the case of an ExpressionSet, \code{exprs(X)} is the data of interest and \code{pData(X)} may contain outcomes and covariates of interest. For most currently implemented tests (exception: tests involving correlation parameters), one hypothesis is tested for each row of the data.}
  \item{W}{A vector or matrix containing non-negative weights to be used in computing the test statistics. If a matrix, \code{W} must be the same dimension as \code{X} with one weight for each value in \code{X}. If a vector, \code{W} may contain one weight for each observation (i.e. column) of \code{X} or one weight for each variable (i.e. row) of \code{X}. In either case, the weights are duplicated appropriately. Weighted F-tests are not available. Default is 'NULL'.}
  \item{Y}{A vector, factor, or \code{Surv} object containing the outcome of interest. This may be class labels (F-tests and two sample t-tests) or a continuous or polycotomous dependent variable (linear regression based t-tests), or survival data (Cox proportional hazards based t-tests). For \code{block.f} and \code{f.twoway} tests, class labels must be ordered by block and within each block ordered by group. If \code{X} is an ExpressionSet, \code{Y} can be a character string referring to the column of \code{pData(X)} to use as outcome. Default is 'NULL'.}
  \item{Z}{A vector, factor, or matrix containing covariate data to be used in the regression (linear and Cox) models. Each variable should be in one column, so that \code{nrow(Z)=ncol(X)}. If \code{X} is an ExpressionSet, \code{Z} can be a character string referring to the column of \code{pData(X)} to use as covariates. The variables \code{Z.incl} and \code{Z.adj} allow one to specify which covariates to use in a particular test without modifying the input \code{Z}. Default is 'NULL'.}
  \item{Z.incl}{The indices of the columns of \code{Z} (i.e. which variables) to include in the model. These can be numbers or column names (if the columns are names). Default is 'NULL'.}
  \item{Z.test}{The index or names of the column of \code{Z} (i.e. which variable) to use to test for association with each row of \code{X} in a linear model. Only used for \code{test="lm.XvsZ"}, where it is necessary to specify which covariate's regression parameter is of interest. Default is 'NULL'.}
  \item{na.rm}{Logical indicating whether to remove observations with an NA. Default is 'TRUE'.}
  \item{test}{Character string specifying the test statistics to use, by default 't.twosamp.unequalvar'. See details (below) for a list of tests.}
  \item{robust}{Logical indicating whether to use the robust version of the chosen test, e.g. Wilcoxon singed rank test for robust one-sample t-test or \code{rlm} instead of \code{lm} in linear models. Default is 'FALSE'.}
  \item{standardize}{Logical indicating whether to use the standardized version of the test statistics (usual t-statistics are standardized). Default is 'TRUE'.}
  \item{alternative}{Character string indicating the alternative hypotheses, by default 'two.sided'. For one-sided tests, use 'less' or 'greater' for null hypotheses of 'greater than or equal' (i.e. alternative is 'less') and 'less than or equal', respectively.}
  \item{psi0}{The hypothesized null value, typically zero (default). Currently, this should be a single value, which is used for all hypotheses.}
  \item{typeone}{Character string indicating which type I error rate to control, by default family-wise error rate ('fwer'). Other options include generalized family-wise error rate ('gfwer'), with parameter \code{k} giving the allowed number of false positives, and tail probability of the proportion of false positives ('tppfp'), with parameter \code{q} giving the allowed proportion of false positives. The false discovery rate ('fdr') can also be controlled.}
  \item{k}{The allowed number of false positives for gFWER control. Default is 0 (FWER).}
  \item{q}{The allowed proportion of false positives for TPPFP control. Default is 0.1.}
  \item{fdr.method}{Character string indicating which FDR controlling method should be used when \code{typeone="fdr"}. The options are "conservative" (default) for the more conservative, general FDR controlling procedure and "restricted" for the method which requires more assumptions.}
  \item{alpha}{The target nominal type I error rate, which may be a vector of error rates. Default is 0.05.}
  \item{smooth.null}{Indicator of whether to use a kernel density estimate for the tail of the null distributon for computing raw pvalues close to zero. Only used if 'rawp' would be zero without smoothing. Default is 'FALSE'.}
  \item{nulldist}{Character string indicating which resampling method to use for estimating the joint test statistics null distribution, by default the non-parametric bootstrap with centering and scaling ('boot.cs').  The old default 'boot' will still compile and will correspond to 'boot.cs'.  Other null distribution options include 'perm', 'boot.ctr', 'boot.qt', and 'ic', corresponding to the permutation distribution, centered-only bootstrap distribution, quantile-transformed bootstrap distribution, and influence curve multivariate normal joint null distribution, respectively.  More details below.}
  \item{B}{The number of bootstrap iterations (i.e. how many resampled data sets), the number of permutations (if \code{nulldist} is 'perm'), or the number of samples from the multivariate normal distribution (if \code{nulldist} is 'ic') Can be reduced to increase the speed of computation, at a cost to precision. Default is 1000.}
  \item{ic.quant.trans}{If \code{nulldist='ic'}, a logical indicating whether or not a marginal quantile transformation using a t-distribution or user-supplied marginal distribution (stored in \code{perm.mat}) should be applied to the multivariate normal null distribution.  Defaults for \code{marg.null} and \code{marg.par} exist, but can also be specified by the user (see below). Default is 'FALSE'.}
  \item{MVN.method}{If \code{nulldist='ic'}, one of 'mvrnorm' or 'Cholesky' designating how correlated normal test statistics are to be generated.  Selecting 'mvrnorm' uses the function of the same name found in the \code{MASS} library, whereas 'Cholesky' relies on a Cholesky decomposition. Default is 'mvrnorm'.}
  \item{penalty}{If \code{nulldist='ic'} and \code{MVN.method='Cholesky'}, the value in \code{penalty} is added to all diagonal elements of the estimated test statistics correlation matrix to ensure that the matrix is positive definite and that internal calls to \code{'chol'} do not return an error.  Default is 1e-6.}
  \item{method}{The multiple testing procedure to use. Options are single-step maxT ('ss.maxT', default), single-step minP ('ss.minP'), step-down maxT ('sd.maxT'), and step-down minP ('sd.minP').}
  \item{get.cr}{Logical indicating whether to compute confidence intervals for the estimates. Not available for F-tests. Default is 'FALSE'.}
  \item{get.cutoff}{Logical indicating whether to compute thresholds for the test statistics. Default is 'FALSE'.}
  \item{get.adjp}{Logical indicating whether to compute adjusted p-values. Default is 'TRUE'.}
  \item{keep.nulldist}{Logical indicating whether to return the computed bootstrap or influence curve null distribution, by default 'TRUE'.  Not available for \code{nulldist}='perm'. Note that this matrix can be quite large.}
  \item{keep.rawdist}{Logical indicating whether to return the computed non-null (raw) bootstrap distribution, by default 'FALSE'.  Not available when using \code{nulldist}='perm' or 'ic'.  Note that this matrix can become quite large.  If one wishes to use subsequent calls to \code{update} or \code{EBupdate} in which one updates choice of bootstrap null distribution, \code{keep.rawdist} must be TRUE.  To save on memory, \code{update} only requires that one of \code{keep.nulldist} or \code{keep.rawdist} be 'TRUE'.}
  \item{seed}{Integer or vector of integers to be used as argument to \code{set.seed} to set the seed for the random number generator for bootstrap resampling. This argument can be used to repeat exactly a test performed with a given seed. If the seed is specified via this argument, the same seed will be returned in the seed slot of the MTP object created. Else a random seed(s) will be generated, used and returned. Vector of integers used to specify seeds for each node in a cluster used to to generate a bootstrap null distribution.}
  \item{cluster}{Integer for number of nodes to create or a cluster object created through the package snow. With \code{cluster=1}, bootstrap is implemented on single node. Supplying a cluster object results in the bootstrap being implemented in parallel on the provided nodes. This option is only available for the bootstrap procedure. With default value of 1, bootstrap is executed on single CPU.}
  \item{type}{Interface system to use for computer cluster. See \code{snow} package for details.}
  \item{dispatch}{The number or percentage of bootstrap iterations to dispatch at a time to each node of the cluster if a computer cluster is used. If dispatch is a percentage, \code{B*dispatch} must be an integer. If dispatch is an integer, then \code{B/dispatch} must be an integer. Default is 5 percent.}
  \item{marg.null}{If \code{nulldist='boot.qt'}, the marginal null distribution to use for quantile transformation.  Can be one of 'normal', 't', 'f' or 'perm'.  Default is 'NULL', in which case the marginal null distribution is selected based on choice of test statistics.  Defaults explained below. If 'perm', the user must supply a vector or matrix of test statistics corresponding to another marginal null distribution, perhaps one created externally by the user, and possibly referring to empirically derived \emph{marginal permutation distributions}, although the statistics could represent any suitable choice of marginal null distribution.}   
  \item{marg.par}{If \code{nulldist='boot.qt'}, the parameters defining the marginal null distribution in \code{marg.null} to be used for quantile transformation.  Default is 'NULL', in which case the values are selected based on choice of test statistics and other available parameters (e.g., sample size, number of groups, etc.).  Defaults explained below.  User can override defaults, in which case a matrix of marginal null distribution parameters can be accepted.  Providing numeric (vector) values will apply the same null distribution defined by the parameter to all hypotheses, while providing a matrix of values allows the user to perform multiple testing using parameters which may vary with each hypothesis, as may be desired in common-quantile minP procedures. In this way, theoretical factors or factors affecting sample size or missingness may be assessed.}
  \item{keep.margpar}{If \code{nulldist='boot.qt'}, a logical indicating whether the (internally created) matrix of marginal null distribution parameters should be returned.  Default is 'TRUE'.}
  \item{ncp}{If \code{nulldist='boot.qt'}, a value for a possible noncentrality parameter to be used during marginal quantile transformation. Default is 'NULL'.}
  \item{perm.mat}{If \code{nulldist='boot.qt'} and \code{marg.null='perm'}, a matrix of user-supplied test statistics from a particular distribution to be used during marginal quantile transformation.  The statistics may represent empirically derived marginal permutation values, may be theoretical values, or may represent a sample from some other suitable choice of marginal null distribution.}
  \item{keep.index}{If \code{nulldist='ic'} and \code{test='t.cor'} or \code{test='z.cor'}, the index returned is a matrix with the indices of the first and second variables considered for pairwise correlations.  If there are p hypotheses, this arguments returns \code{t(combn(p,2))}.  For all other choices of test statistic, the index is not returned, as they correspond to the original order of the hypotheses in \code{X}.}  
  \item{keep.label}{Default is 'FALSE'.  A logical indicating whether or not the label in \code{Y} should be returned as a slot in the resulting MTP object.  Typically not necessary, although useful if one is using \code{update} and wants to use marginal null distribution defaults with \code{nulldist='boot.qt'} (e.g., with F-tests).}
}

\details{
A multiple testing procedure (MTP) is defined by choices of test statistics, type I error rate, null distribution and method for error rate control. Each component is described here. For two-sample t-tests, the group with the smaller-valued label is substracted from the group with the larger-valued label.  That is, differences in means are calculated as "mean of group 2 - mean of group 1" or "mean of group B - mean of group A". For paired t-tests, the arrangement of group indices does not matter, as long as the columns are arranged in the same corresponding order between groups. For example, if group 1 is coded as 0, and group 2 is coded as 1, for 3 pairs of data, it does not matter if the label \code{Y} is coded as "0,0,0,1,1,1", "1,1,1,0,0,0" "0,1,0,1,0,1" or "1,0,1,0,1,0", the paired differences between groups will be calculated as "group 2 - group 1". See references for more detail. 

Test statistics are determined by the values of \code{test}: 
\describe{
\item{t.onesamp:}{one-sample t-statistic for tests of means;}
\item{t.twosamp.equalvar:}{equal variance two-sample t-statistic for tests of differences in means (two-sample t-statistic);}
\item{t.twosamp.unequalvar:}{unequal variance two-sample t-statistic for tests of differences in means (two-sample Welch t-statistic);}
\item{t.pair:}{two-sample paired t-statistic for tests of differences in means;}
\item{f:}{multi-sample F-statistic for tests of equality of population means (assumes constant variance across groups, but not normality); }
\item{f.block:}{multi-sample F-statistic for tests of equality of population means in a block design (assumes constant variance across groups, but not normality). This test is not available with the bootstrap null distribution;}
\item{f.twoway:}{multi-sample F-statistic for tests of equality of population means in a block design (assumes constant variance across groups, but not normality). Differs from \code{f.block} in requiring multiple observations per group*block combintation. This test uses the means of each group*block combination as response variable and test for group main effects assuming a randomized block design;}
\item{lm.XvsZ:}{t-statistic for tests of regression coefficients for variable \code{Z.test} in linear models, each with a row of X as outcome, possibly adjusted by covariates \code{Z.incl} from the matrix \code{Z} (in the case of no covariates, one recovers the one-sample t-statistic, \code{t.onesamp});}
\item{lm.YvsXZ:}{t-statistic for tests of regression coefficients in linear models, with outcome Y and each row of X as covariate of interest, with possibly other covariates \code{Z.incl} from the matrix \code{Z};}
\item{coxph.YvsXZ:}{t-statistic for tests of regression coefficients in Cox proportional hazards survival models, with outcome Y and each row of X as covariate of interest, with possibly other covariates \code{Z.incl} from the matrix \code{Z}.}
\item{t.cor}{t-statistics for tests of pairwise correlation parameters for all variables in X.  Note that the number of hypotheses can become quite large very fast.  This test is only available with the influence curve null distribution.}
\item{z.cor}{Fisher's z-statistics for tests of pairwise correlation parameters for all variables in X.  Note that the number of hypotheses can become quite large very fast.  This test is only available with the influence curve null distribution.}
}

When \code{robust=TRUE}, non-parametric versions of each test are performed. For the linear models, this means \code{rlm} is used instead of \code{lm}. There is not currently a robust version of \code{test=coxph.YvsXZ}. For the t- and F-tests, data values are simply replaced by their ranks. This is equivalent to performing the following familiar named rank-based tests. The conversion after each test is the formula to convert from the MTP test to the statistic reported by the listed R function (where num is the numerator of the MTP test statistics, n is total sample size, nk is group k sample size, K is total number of groups or treatments, and rk are the ranks in group k).
\describe{
\item{t.onesamp or t.pair:}{Wilcoxon signed rank, \code{wilcox.test} with \code{y=NULL} or \code{paired=TRUE}, \cr
conversion: num/n}
\item{t.twosamp.equalvar:}{Wilcoxon rank sum or Mann-Whitney, \code{wilcox.test}, \cr
conversion: n2*(num+mean(r1)) - n2*(n2+1)/2}
\item{f:}{Kruskal-Wallis rank sum, \code{kruskal.test}, \cr
conversion: num*12/(n*(n-1))}
\item{f.block:}{Friedman rank sum, \code{friedman.test}, \cr
conversion: num*12/(K*(K+1))}
\item{f.twoway:}{Friedman rank sum, \code{friedman.test}, \cr
conversion: num*12/(K*(K+1))}
}

The implemented MTPs are based on control of the family-wise error rate, defined as the probability of any false positives. Let Vn denote the (unobserved) number of false positives. Then, control of FWER at level alpha means that Pr(Vn>0)<=alpha. The set of rejected hypotheses under a FWER controlling procedure can be augmented to increase the number of rejections, while controlling other error rates. The generalized family-wise error rate is defined as Pr(Vn>k)<=alpha, and it is clear that one can simply take an FWER controlling procedure, reject k more hypotheses and have control of gFWER at level alpha. The tail probability of the proportion of false positives depends on both the number of false postives (Vn) and the number of rejections (Rn). Control of TPPFP at level alpha means Pr(Vn/Rn>q)<=alpha, for some proportion q. Control of the false discovery rate refers to the expected proportion of false positives (rather than a tail probability). Control of FDR at level alpha means E(Vn/Rn)<=alpha.

In practice, one must choose a method for estimating the test statistics null distribution. We have implemented several versions of an ordinary non-parametric bootstrap estimator and a permutation estimator (which makes sense in certain settings, see references). The non-parametric bootstrap estimator (default) provides asymptotic control of the type I error rate for any data generating distribution, whereas the permutation estimator requires the subset pivotality assumption. One draw back of both methods is the discreteness of the estimated null distribution when the sample size is small. Furthermore, when the sample size is small enough, it is possible that ties will lead to a very small variance estimate. Using \code{standardize=FALSE} allows one to avoid these unusually small test statistic denominators. Parametric bootstrap estimators are another option (not yet implemented).  For asymptotically linear estimators, such as those commonly probed using t-statistics, another choice of null distribution is provided when sampling from a multivariate normal distribution with mean zero and correlation matrix derived from the vector influence function.  Sampling from a multivariate normal may alleviate the discreteness of the bootstrap and permutation distributions, although accuracy in estimation of the test statistics correlation matrix will be of course also affected by sample size.  

For the nonparametric bootstrap distribution with marginal null quantile transformation, the following defaults for \code{marg.null} and \code{marg.par} are available based on choice of test statistics, sample size 'n', and various other parameters:
\describe{
\item{t.onesamp:}{t-distribution with df=n-1;}
\item{t.twosamp.equalvar:}{t-distribution with df=n-2;}
\item{t.twosamp.unequalvar:}{N(0,1);}
\item{t.pair:}{t-distribution with df=n-1, where n is the number of unique samples, i.e., the number of observed differences between paired samples;}
\item{f:}{F-distribution with df1=k-1, df2=n-k, for k groups;}
\item{f.block:}{NA. Only available with permutation distribution;}
\item{f.twoway:}{F-distribution with df1=k-1,df2=n-k*l, for k groups and l blocks;}
\item{lm.XvsZ:}{N(0,1);}
\item{lm.YvsXZ:}{N(0,1);}
\item{coxph.YvsXZ:}{N(0,1);}
\item{t.cor}{t-distribution with df=n-2;}
\item{z.cor}{N(0,1).}
}

The above defaults, however, can be overridden by manually setting values of \code{marg.null} and \code{marg.par}.  In the case of \code{nulldist='ic'}, and \code{ic.quant.trans=TRUE}, the defaults are the same as above except that 'lm.XvsZ' and 'lm.YvsXZ' are replaced with t-distributions with df=n-p.
  
Given observed test statistics, a type I error rate (with nominal level), and a test statistics null distribution, MTPs provide adjusted p-values, cutoffs for test statistics, and possibly confidence regions for estimates. Four methods are implemented, based on minima of p-values and maxima of test statistics. Only the step down methods are currently available with the permutation null distribution.

Computation times using a bootstrap null distribution are slower when weights are used for one and two-sample tests. Computation times when using a bootstrap null distribution also are slower for the tests \code{lmXvsZ}, \code{lmYvsXZ}, \code{coxph.YvsXZ}.

To execute the bootstrap on a computer cluster, a cluster object generated with \code{makeCluster} in the package \code{snow} may be used as the argument for cluster. Alternatively, the number of nodes to use in the computer cluster can be used as the argument to cluster. In this case, \code{type} must be specified and a cluster will be created. In both cases, \code{Biobase} and \code{multtest} will be loaded onto each cluster node if these libraries are located in a directory in the standard search path. If these libraries are in a non-standard location, it is necessary to first create the cluster, load \code{Biobase} and \code{multtest} on each node and then to use the cluster object as the argument to cluster. See documentation for \code{snow} package for additional information on creating and using a cluster.

Finally, note that the old argument \code{csnull} is now DEPRECATED as of \code{multtest} v. 2.0.0 given the expanded null distribution options described above.  Previously, this argument was an indicator of whether the bootstrap estimated test statistics distribution should be centered and scaled (to produce a null distribution) or not. If \code{csnull=FALSE}, the (raw) non-null bootstrap estimated test statistics distribution was returned.  If the non-null bootstrap distribution should be returned, this object is now stored in the 'rawdist' slot when \code{keep.rawdist=TRUE} in the original \code{MTP} function call.  
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
   \item{\code{reject}}{Object of class \code{'matrix'}, rejection indicators (TRUE for a rejected null hypothesis), for each value of the nominal Type I error rate \code{alpha}.}
    \item{\code{rawdist}}{The numeric matrix for the estimated nonparametric non-null test statistics distribution (returned only if \code{keep.rawdist=TRUE} and if \code{nulldist} is one of 'boot.ctr', 'boot.cs', or 'boot.qt'). This slot must not be empty if one wishes to call \code{update} to change choice of bootstrap-based null distribution.}
    \item{\code{nulldist}}{The numeric matrix for the estimated test statistics null distribution (returned only if \code{keep.nulldist=TRUE}); option not currently available for permutation null distribution, i.e.,  \code{nulldist='perm'}). By default (i.e., for \code{nulldist='boot.cs'}), the entries of \code{nulldist} are the null value shifted and scaled bootstrap test statistics, with one null test statistic value for each hypothesis (rows) and bootstrap iteration (columns).}
    \item{\code{nulldist.type}}{Character value describing which choice of null distribution was used to generate the MTP results.  Takes on one of the values of the original \code{nulldist} argument in the call to MTP, i.e., 'boot.cs', 'boot.ctr', 'boot.qt', 'ic', or 'perm'.}
   \item{\code{marg.null}}{If \code{nulldist='boot.qt'}, a character value returning which choice of marginal null distribution was used by the MTP.  Can be used to check default values or to ensure manual settings were correctly applied.}
   \item{\code{marg.par}}{If \code{nulldist='boot.qt'}, a numeric matrix returning the parameters of the marginal null distribution(s) used by the MTP.  Can be used to check default values or to ensure manual settings were correctly applied.}
   \item{\code{call}}{Object of class \code{call}, the call to the MTP function.}
    \item{\code{seed}}{An integer or vector for specifying the state of the random number generator used to create the resampled datasets. The seed can be reused for reproducibility in a repeat call to \code{MTP}. This argument is currently used only for the bootstrap null distribution (i.e., for \code{nulldist="boot.xx"}). See \code{?set.seed} for details.}
}

\references{
M.J. van der Laan, S. Dudoit, K.S. Pollard (2004), Augmentation Procedures for Control of the Generalized Family-Wise Error Rate and Tail Probabilities for the Proportion of False Positives, Statistical Applications in Genetics and Molecular Biology, 3(1). 
\url{http://www.bepress.com/sagmb/vol3/iss1/art15/}

M.J. van der Laan, S. Dudoit, K.S. Pollard (2004), Multiple Testing. Part II. Step-Down Procedures for Control of the Family-Wise Error Rate, Statistical Applications in Genetics and Molecular Biology, 3(1).
\url{http://www.bepress.com/sagmb/vol3/iss1/art14/}

S. Dudoit, M.J. van der Laan, K.S. Pollard (2004), Multiple Testing. Part I. Single-Step Procedures for Control of General Type I Error Rates, Statistical Applications in Genetics and Molecular Biology, 3(1).
\url{http://www.bepress.com/sagmb/vol3/iss1/art13/}

K.S. Pollard and Mark J. van der Laan, "Resampling-based Multiple Testing: Asymptotic Control of Type I Error and Applications to Gene Expression Data" (June 24, 2003). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 121.
\url{http://www.bepress.com/ucbbiostat/paper121}

M.J. van der Laan and A.E. Hubbard (2006), Quantile-function Based Null Distributions in Resampling Based Multiple Testing, Statistical Applications in Genetics and Molecular Biology, 5(1).
\url{http://www.bepress.com/sagmb/vol5/iss1/art14/}

S. Dudoit and M.J. van der Laan.  Multiple Testing Procedures and Applications to Genomics.  Springer Series in Statistics. Springer, New York, 2008.

}

\author{Katherine S. Pollard and Houston N. Gilbert with design contributions from Sandra Taylor, Sandrine Dudoit and Mark J. van der Laan.}

\note{Thank you to Peter Dimitrov for suggestions about the code.}

\seealso{\code{\link{EBMTP}}, \code{\link{MTP-class}}, \code{\link{MTP-methods}}, \code{\link{mt.minP}}, \code{\link{mt.maxT}}, \code{\link{ss.maxT}}, \code{\link{fwer2gfwer}}}

\examples{
 
#data 
set.seed(99)
data<-matrix(rnorm(90),nr=9)
group<-c(rep(1,5),rep(0,5))

#fwer control with centered and scaled bootstrap null distribution 
#(B=100 for speed)
m1<-MTP(X=data,Y=group,alternative="less",B=100,method="sd.minP")
print(m1)
summary(m1)
par(mfrow=c(2,2))
plot(m1,top=9)

#fwer control with quantile transformed bootstrap null distribution
#default settings = N(0,1) marginal null distribution
m2<-MTP(X=data,Y=group,alternative="less",B=100,method="sd.minP",
	nulldist="boot.qt",keep.rawdist=TRUE)

#fwer control with quantile transformed bootstrap null distribution
#marginal null distribution and df parameters manually set, 
#first all equal, then varying with hypothesis
m3<-update(m2,marg.null="t",marg.par=10)
mps<-matrix(c(rep(9,5),rep(10,5)),nr=10,nc=1)
m4<-update(m2,marg.null="t",marg.par=mps)

m1@nulldist.type
m2@nulldist.type
m2@marg.null
m2@marg.par
m3@nulldist.type
m3@marg.null
m3@marg.par
m4@nulldist.type
m4@marg.null
m4@marg.par

}

\keyword{htest}


