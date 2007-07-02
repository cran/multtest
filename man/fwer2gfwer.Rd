\name{fwer2gfwer}

\alias{fwer2gfwer}
\alias{fwer2tppfp}
\alias{fwer2fdr}

\title{Function to compute augmentation MTP adjusted p-values}

\description{Augmentation multiple testing procedures (AMTPs) to control the generalized family-wise error rate (gFWER), the tail probability of the proportion of false positives (TPPFP), and false discovery rate (FDR) based on any initial procudeure controlling the family-wise error rate (FWER). AMTPs are obtained by adding suitably chosen null hypotheses to the set of null hypotheses already rejected by an initial FWER-controlling MTP. A function for control of FDR given any TPPFP controlling procedure is also provided.
}

\usage{
fwer2gfwer(adjp, k = 0)

fwer2tppfp(adjp, q = 0.05)

fwer2fdr(adjp, method = "both", alpha = 0.05)

}


\arguments{
  \item{adjp}{Numeric vector of adjusted p-values from any FWER-controlling procedure.}
  \item{k}{Maximum number of false positives.}
  \item{q}{Maximum proportion of false positives.}
  \item{method}{Character string indicating which FDR controlling method should be used. The options are "conservative" for a conservative, general method, "restricted" for a less conservative, but restricted method, or "both" (default) for both.}
  \item{alpha}{Nominal level for an FDR controlling procedure (can be a vector of levels).}
}

\details{
The gFWER and TPPFP functions control Type I error rates defined as tail probabilities for functions g(Vn,Rn) of the numbers of Type I errors (Vn) and rejected hypotheses (Rn). The gFWER and TPPFP correspond to the special cases g(Vn,Rn)=Vn (number of false positives) and g(Vn,Rn)=Vn/Rn (proportion of false positives among the rejected hypotheses), respectively. 

Adjusted p-values for an AMTP are simply shifted versions of the adjusted p-values of the original FWER-controlling MTP. For control of gFWER (Pr(Vn>k)), for example, the first \code{k} adjusted p-values are set to zero and the remaining p-values are the adjusted p-values of the FWER-controlling MTP shifted by k. One can therefore build on the large pool of available FWER-controlling procedures, such as the single-step and step-down maxT and minP procedures.

Given a FWER-controlling MTP, the FDR can be conservatively controlled at level \code{alpha} by considering the corresponding TPPFP AMTP with \code{q=alpha/2} at level \code{alpha/2}, so that Pr(Vn/Rn>alpha/2)<=alpha/2. A less conservative procedure (\code{general=FALSE}) is obtained by using an AMTP controlling the TPPFP with \code{q=1-sqrt(1-alpha)} at level \code{1-sqrt(1-alpha)}, so that Pr(Vn/Rn>1-sqrt(1-alpha))<=1-sqrt(1-alpha). The first, more general method can be used with any procedure that asymptotically controls FWER. The second, less conservative method requires the following additional assumptions: (i) the true alternatives are asymptotically always rejected by the FWER-controlling procedure, (ii) the limit of the FWER exists, and (iii) the FWER-controlling procedure provides exact asymptotic control. See \url{http://www.bepress.com/sagmb/vol3/iss1/art15/} for more details. The method implemented in \code{fwer2fdr} for computing rejections simply uses the TPPFP AMTP \code{fwer2tppfp} with \code{q=alpha/2} (or 1-sqrt(1-alpha)) and rejects each hypothesis for which the TPPFP adjusted p-value is less than or equal to alpha/2 (or 1-sqrt(1-alpha)). The adjusted p-values are based directly on the FWER adjusted p-values, so that very occasionally a hypothesis will have the indicator that it is rejected in the matrix of rejections, but the adjusted p-value will be slightly greater than the nominal level. The opposite might also occur occasionally.
}

\value{
For \code{fwer2gfwer} and \code{fwer2tppfp}, a numeric vector of AMTP adjusted p-values. For \code{fwer2fdr}, a list with two components: (i) a numeric vector (or a \code{length(adjp)} by 2 matrix if \code{method="both"}) of adjusted p-values for each hypothesis, (ii) a \code{length(adjp)} by \code{length(alpha)} matrix (or \code{length(adjp)} by \code{length(alpha)} by 2 array if \code{method="both"}) of indicators of whether each hypothesis is rejected at each value of the argument \code{alpha}.
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

\author{Katherine S. Pollard with design contributions from Sandrine Dudoit and Mark J. van der Laan.}

\seealso{\code{\link{MTP}}, \code{\link{MTP-class}}, \code{\link{MTP-methods}}, \code{\link{mt.minP}}, \code{\link{mt.maxT}}}

\examples{

data<-matrix(rnorm(200),nr=20)
group<-c(rep(0,5),rep(1,5))
fwer.mtp<-MTP(X=data,Y=group)
fwer.adjp<-fwer.mtp@adjp
gfwer.adjp<-fwer2gfwer(adjp=fwer.adjp,k=c(1,5,10))
compare.gfwer<-cbind(fwer.adjp,gfwer.adjp)
mt.plot(adjp=compare.gfwer,teststat=fwer.mtp@statistic,proc=c("gFWER(0)","gFWER(1)","gFWER(5)","gFWER(10)"),col=1:4,lty=1:4)
title("Comparison of Single-step MaxT gFWER Controlling Methods")

}

\keyword{htest}
\keyword{internal}



