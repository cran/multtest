\name{multtest-internal}
\alias{.mt.BLIM}
\alias{.mt.RandSeed}
\alias{.mt.naNUM}
\alias{mt.number2na}
\alias{mt.na2number}
\alias{mt.getmaxB}
\alias{mt.transformL}
\alias{mt.transformX}
\alias{mt.transformV}
\alias{mt.checkothers}
\alias{mt.checkX}
\alias{mt.checkV}
\alias{mt.checkclasslabel}
\alias{mt.niceres}
\alias{mt.legend}
\title{Internal multtest functions and variables}
\description{
  Internal multtest functions and variables
}
\usage{
.mt.BLIM
.mt.RandSeed
.mt.naNUM
mt.number2na(x,na)
mt.na2number(x,na)
mt.getmaxB(classlabel,test,B, verbose)
mt.transformL(classlabel,test)
mt.transformV(V,classlabel,test,na,nonpara)
mt.transformX(X,classlabel,test,na,nonpara)
mt.checkothers(side="abs",fixed.seed.sampling="y", B=10000,
na=.mt.naNUM, nonpara="n")
mt.checkX(X,classlabel,test)
mt.checkV(V,classlabel,test)
mt.checkclasslabel(classlabel,test)
mt.niceres<-function(res,X,index)
mt.legend(x, y = NULL, legend, fill = NULL, col = "black", lty, 
    lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"), 
    pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, xjust = 0, 
    yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), 
    text.width = NULL, text.col = par("col"), merge = do.lines && 
        has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE,...)
}
\details{
  These are not to be called directly by the user.
}
\keyword{internal}
