\name{MTP-methods}
\docType{methods}
\alias{MTP-methods}
\alias{EBMTP-methods}

\alias{[-methods}
\alias{[,MTP-method}
\alias{[,EBMTP-method}

\alias{as.list-methods}
\alias{as.list,MTP-method}
\alias{as.list,EBMTP-method}
\alias{as.list}

\alias{plot-methods}
\alias{plot,MTP-method}
\alias{plot,EBMTP-method}
\alias{plot}

\alias{print-methods}
\alias{print,MTP-method}
\alias{print,EBMTP-method}
\alias{print.MTP}

\alias{summary-methods}
\alias{summary,MTP-method} 
\alias{summary,EBMTP-method} 
\alias{summary}

\alias{update-methods}
\alias{update,MTP-method}
\alias{update}

\alias{EBupdate-methods}
\alias{EBupdate,EBMTP-method}
\alias{EBupdate}

\alias{mtp2ebmtp-methods}
\alias{mtp2ebmtp,MTP-method}
\alias{mtp2ebmtp}

\alias{ebmtp2mtp-methods}
\alias{ebmtp2mtp,EBMTP-method}
\alias{ebmtp2mtp}


\title{Methods for MTP and EBMTP objects in Package `multtest'}

\description{Summary, printing, plotting, subsetting, updating, as.list and class conversion methods were defined for the \code{MTP} and \code{EBMTP} classes. These methods provide visual and numeric summaries of the results of a multiple testing procedure (MTP) and allow one to perform some basic manipulations of  objects class \code{MTP} or \code{EBMTP}.  \cr

Several of the methods with the same name will work on objects of their respective class.  One exception to this rule is the difference between \code{update} and \code{EBupdate} (described below).  Because of the differences in the testing procedures, separately named methods were chosen to clearly delineate which method was being applied to which type of object.}

\section{Methods}{
\describe{

   \item{[}{: Subsetting method for \code{MTP} and \code{EBMTP} classes, which operates selectively on each slot of an \code{MTP} or \code{EBMTP} instance to retain only the data related to the specified hypotheses.}

    \item{as.list}{: Converts an object of class \code{MTP} or \code{EBMTP} to an object of class \code{list}, with an entry for each slot.}
 
   \item{plot}{: plot methods for \code{MTP} and \code{EBMTP} classes, produces the following graphical summaries of the results of a MTP. The type of display may be specified via the \code{which} argument. \cr

1. Scatterplot of number of rejected hypotheses vs. nominal Type I error rate. \cr

2. Plot of ordered adjusted p-values; can be viewed as a plot of Type I error rate vs. number of rejected hypotheses. \cr

3. Scatterplot of adjusted p-values vs. test statistics (also known as volcano plot). \cr

4. Plot of unordered adjusted p-values. \cr

Only for objects of class \code{MTP}: \cr

5. Plot of confidence regions for user-specified parameters, by default the 10 parameters corresponding to the smallest adjusted p-values  (argument \code{top}). \cr

6. Plot of test statistics and corresponding cut-offs (for each value of \code{alpha}) for user-specified hypotheses, by default the 10 hypotheses corresponding to the smallest adjusted p-values (argument \code{top}). \cr

Plots (5) and (6) are not available for objects of class \code{EBMTP} because the function \code{EBMTP} returns only adjusted p-values and not confidence regions of cut-offs.  The argument \code{logscale} (by default equal to FALSE) allows one to use the negative decimal logarithms of the adjusted p-values in the second, third, and fourth graphical displays. The arguments \code{caption} and \code{sub.caption} allow one to change the titles and subtitles for each of the plots (default subtitle is the MTP function call). Note that some of these plots are implemented in the older function \code{mt.plot}.}

    \item{print}{: print method for \code{MTP} and \code{EBMTP} classes, returns a description of an object of either class, including sample size, number of tested hypotheses, type of test performed (value of argument \code{test}), Type I error rate (value of argument \code{typeone}), nominal level of the test (value of argument \code{alpha}), name of the MTP (value of argument \code{method}), call to the function \code{MTP} or \code{EBMTP}. 

In addition, this method produces a table with the class, mode, length, and dimension of each slot of the \code{MTP} or \code{EBMTP} instance. 
}

    \item{summary}{: summary method for \code{MTP} and \code{EBMTP} classes, provides numerical summaries of the results of a MTP and returns a list with the following three components. \cr

1. rejections: A data.frame with the number(s) of rejected hypotheses for the nominal Type I error rate(s) specified by the \code{alpha} argument of the function \code{MTP} or \code{EBMTP}. (For objects of class \code{MTP}, NULL values are returned if all three arguments \code{get.cr}, \code{get.cutoff}, and \code{get.adjp} are FALSE). \cr

2. index: A numeric vector of indices for ordering the hypotheses according to first \code{adjp}, then \code{rawp}, and finally the absolute value of \code{statistic} (not printed in the summary). \cr 

3. summaries: When applicable (i.e., when the corresponding quantities are returned by \code{MTP} or \code{EBMTP}), a table with six number summaries of the distributions of the adjusted p-values, unadjusted p-values, test statistics, and parameter estimates.}
    \item{update}{: update methods for \code{MTP} class, respectively, provides a mechanism to re-run the MTP with different choices of the following arguments - nulldist, alternative, typeone, k, q, fdr.method, alpha, smooth.null, method, get.cr, get.cutoff, get.adjp, keep.nulldist, keep.rawdist, keep.margpar. When evaluate is 'TRUE', a new object of class MTP is returned. Else, the updated call is returned. The \code{MTP} object passed to the update method must have either a  non-empty \code{rawdist} slot or a non-empty \code{nulldist} slot (i.e., must have been called with either 'keep.rawdist=TRUE' or 'keep.nulldist=TRUE').\cr}

    \item{EBupdate}{: update method for \code{EBMTP} class, provides a mechanism to re-run the MTP with different choices of the following arguments - nulldist, alternative, typeone, k, q, alpha, smooth.null, bw, kernel, prior, keep.nulldist, keep.rawdist, keep.falsepos, keep.truepos, keep.errormat, keep.margpar. When evaluate is 'TRUE', a new object of class EBMTP is returned. Else, the updated call is returned. The \code{EBMTP} object passed to the update method must have either a  non-empty \code{rawdist} slot or a non-empty \code{nulldist} slot (i.e., must have been called with either 'keep.rawdist=TRUE' or 'keep.nulldist=TRUE').\cr 

Additionally, when calling \code{EBupdate} for any Type I error rate other than FWER, the \code{typeone} argument must be specified (even if the original object did not control FWER). For example, 
\code{typeone="fdr"}, would always have to be specified, even if the original object also controlled the FDR. In other words, for all function arguments, it is safest to always assume that you 
are updating from the \code{EBMTP} default function settings, regardless of the original call to the \code{EBMTP} function. Currently, the main advantage of the \code{EBupdate} method is that it prevents the need for repeated estimation of the test statistics null distribution. \cr 

To save on memory, if one knows ahead of time that one will want to compare different choices of bootstrap-based null distribution, then it is both necessary and sufficient to specify 'keep.rawdist=TRUE', as there is no other means of moving between null distributions other than through the non-transformed non-parametric bootstrap distribution.  In this case, 'keep.nulldist=FALSE' may be used.  Specifically, if an object of class \code{MTP} or \code{EBMTP} contains a non-empty \code{rawdist} slot and an empty \code{nulldist} slot, then a new null distribution will be generated according to the values of the \code{nulldist=} argument in the original call to \code{(EB)MTP} or any additional specifications in the call to \code{(EB)update}. On the other hand, if one knows that one wishes to only update an \code{(EB)MTP} object in ways which do not involve choice of bootstrap null distribution, then 'keep.nulldist=TRUE' will suffice and 'keep.rawdist' can be set to \code{FALSE} (default settings).  The original null distribution object will then be used for all subsequent calls to \code{update}. \cr

N.B.: Note that \code{keep.rawdist=TRUE} is only available for the bootstrap-based resampling methods.  The non-null distribution does not exist for the permutation or influence curve multivariate normal null distributions.
}
    \item{mtp2ebmtp}{: coersion method for converting objects of class \code{MTP} to objects of class \code{EBMTP}.  Slots common to both objects are taken from the object of class \code{MTP} and used to create a new object of class \code{EBMTP}.  Once an object of class \code{EBMTP} is created, one may use the method \code{EBupdate} to perform resampling-based empirical Bayes multiple testing without the need for repeated resampling.}
    \item{ebmtp2mtp}{: coersion method for converting objects of class \code{EBMTP} to objects of class \code{MTP}.  Slots common to both objects are taken from the object of class \code{EBMTP} and used to create a new object of class \code{MTP}.  Once an object of class \code{MTP} is created, one may use the method \code{update} to perform resampling-based multiple testing (as would have been done with calls to \code{MTP}) without the need for repeated resampling.}
}
}


\author{Katherine S. Pollard and Houston N. Gilbert with design contributions from Sandrine Dudoit and Mark J. van der Laan.}

\keyword{methods}




