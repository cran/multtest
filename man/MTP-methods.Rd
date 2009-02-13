\name{MTP-methods}
\docType{methods}
\alias{MTP-methods}

\alias{[-methods}
\alias{[,MTP-method}

\alias{as.list-methods}
\alias{as.list,MTP-method}
\alias{as.list}

\alias{plot-methods}
\alias{plot,MTP-method}
\alias{plot}

\alias{print-methods}
\alias{print,MTP-method}
\alias{print.MTP}

\alias{summary-methods}
\alias{summary,MTP-method} 
\alias{summary}

\alias{update-methods}
\alias{update,MTP-method}
\alias{update}

\title{Methods for MTP objects in Package `multtest'}

\description{Summary, printing, plotting, subsetting, updating and as.list methods were defined for the \code{MTP} class. These methods provide visual and numeric summaries of the results of a multiple testing procedure (MTP) and allow one to perform some basic manipulations of an object of class \code{MTP}.}

\section{Methods}{
\describe{

   \item{[}{: Subsetting method for \code{MTP} class, which operates selectively on each slot of an \code{MTP} instance to retain only the data related to the specified hypotheses.}

    \item{as.list}{: Converts an object of class \code{MTP} to an object of class \code{list}, with an entry for each slot.}
 
   \item{plot}{: plot methods for \code{MTP} class, produces the following graphical summaries of the results of a MTP. The type of display may be specified via the \code{which} argument. \cr
1. Scatterplot of number of rejected hypotheses vs. nominal Type I error rate. \cr

2. Plot of ordered adjusted p-values; can be viewed as a plot of Type I error rate vs. number of rejected hypotheses. \cr

3. Scatterplot of adjusted p-values vs. test statistics (also known as volcano plot). \cr

4. Plot of unordered adjusted p-values. \cr

5. Plot of confidence regions for user-specified parameters, by default the 10 parameters corresponding to the smallest adjusted p-values  (argument \code{top}). \cr

6. Plot of test statistics and corresponding cut-offs (for each value of \code{alpha}) for user-specified hypotheses, by default the 10 hypotheses corresponding to the smallest adjusted p-values (argument \code{top}). \cr

The argument \code{logscale} (by default equal to FALSE) allows one to use the negative decimal logarithms of the adjusted p-values in the second, third, and fourth graphical displays. The arguments \code{caption} and \code{sub.caption} allow one to change the titles and subtitles for each of the plots (default subtitle is the MTP function call). Note that some of these plots are implemented in the older function \code{mt.plot}.}

    \item{print}{: print method for \code{MTP} class, returns a description of an object of class \code{MTP}, including sample size, number of tested hypotheses, type of test performed (value of argument \code{test}), Type I error rate (value of argument \code{typeone}), nominal level of the test (value of argument \code{alpha}), name of the MTP  (value of argument \code{method}), call to the function \code{MTP}. 

In addition, this method produces a table with the class, mode, length, and dimension of each slot of the \code{MTP} instance. 
}

    \item{summary}{: summary method for \code{MTP} class, provides numerical summaries of the results of a MTP and returns a list with the following three components. \cr

1. rejections: A data.frame with the number(s) of rejected hypotheses for the nominal Type I error rate(s) specified by the \code{alpha} argument of the function \code{MTP}. (NULL values are returned if all three arguments \code{get.cr}, \code{get.cutoff}, and \code{get.adjp} are FALSE). \cr

2. index: A numeric vector of indices for ordering the hypotheses according to first \code{adjp}, then \code{rawp}, and finally the absolute value of \code{statistic} (not printed in the summary). \cr 

3. summaries: When applicable (i.e., when the corresponding quantities are returned by \code{MTP}), a table with six number summaries of the distributions of the adjusted p-values, unadjusted p-values, test statistics, and parameter estimates.}
    \item{update}{: update method for \code{MTP} class, provides a mechanism to re-run the MTP with different choices of the following arguments - alternative, typeone, k, q, fdr.method, alpha, smooth.null, method, get.cr, get.cutoff, get.adjp, keep.nulldist. When evaluate is 'TRUE', a new object of class MTP is returned. Else, the updated call is returned. The \code{MTP} object passed to the update method must have a non-empty nulldist slot (ie: must have been called with 'keep.nulldist=TRUE').}
} 
}

\author{Katherine S. Pollard with design contributions from Sandrine Dudoit and Mark J. van der Laan.}

\keyword{methods}




