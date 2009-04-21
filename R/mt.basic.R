mt.rawp2adjp<-function(rawp,
      proc=c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY","ABH",
      "TSBH"), alpha=0.05)
{

  m<-length(rawp)
  n<-length(proc)
  a<-length(alpha)
  index<-order(rawp)
  h0.ABH<-NULL
  h0.TSBH <- NULL
  spval<-rawp[index]

  adjp<-matrix(0,m,n+1)
  dimnames(adjp)<-list(NULL,c("rawp",proc))
  adjp[,1]<-spval

  if(is.element("TSBH",proc))
  {
    #N.B.: This method performed first in order to handle a potential $adjp
    #dimension change in the case that length(alpha)>1.
    #Could also be possibly done using more append() functions, should more 
    #alpha-dependent procedures be developed/included later.
    TS.spot <- which(proc=="TSBH")
    TSBHs<-paste("TSBH",alpha,sep="_")
    newprocs<-append(proc,TSBHs,after=TS.spot)
    newprocs<-newprocs[newprocs!="TSBH"]
    adjp<-matrix(0,m,n+a)
    dimnames(adjp)<-list(NULL,c("rawp",newprocs))
    adjp[,1]<-spval

    # Apply first-pass BH.  
    tmp<-spval
    for(i in (m-1):1){
      tmp[i]<-min(tmp[i+1],min((m/i)*spval[i],1,na.rm=TRUE),na.rm=TRUE)
      if(is.na(spval[i])) tmp[i]<-NA
    }
    # Now use first-pass results to estimate h_0, the number of true nulls.
    # These results depend on the nominal testing level, alpha.
    h0.TSBH <- rep(0,length(alpha))
    names(h0.TSBH) <- paste("h0.TSBH",alpha,sep="_")
    for(i in 1:length(alpha)){
      h0.TSBH[i] <- m - sum(tmp < alpha[i]/(1+alpha[i]),na.rm=TRUE)
      adjp[,TS.spot+i]<-tmp*h0.TSBH[i]/m
    }
  }

  if(is.element("Bonferroni",proc))
  {
    tmp<-m*spval
    tmp[tmp>1]<-1
    adjp[,"Bonferroni"]<-tmp
  }

  if(is.element("Holm",proc))
  {
    tmp<-spval
    tmp[1]<-min(m*spval[1],1)
    for(i in 2:m)
      tmp[i]<-max(tmp[i-1],min((m-i+1)*spval[i],1))
    adjp[,"Holm"]<-tmp
  }

  if(is.element("Hochberg",proc))
  {
    tmp<-spval
    for(i in (m-1):1){
      tmp[i]<-min(tmp[i+1],min((m-i+1)*spval[i],1,na.rm=TRUE),na.rm=TRUE)
      if(is.na(spval[i])) tmp[i]<-NA
    }
    adjp[,"Hochberg"]<-tmp
  }

  if(is.element("SidakSS",proc))
    adjp[,"SidakSS"]<-1-(1-spval)^m

  if(is.element("SidakSD",proc))
  {
    tmp<-spval
    tmp[1]<-1-(1-spval[1])^m
    for(i in 2:m)
      tmp[i]<-max(tmp[i-1],1-(1-spval[i])^(m-i+1))
    adjp[,"SidakSD"]<-tmp
  }

  if(is.element("BH",proc))
  {
    tmp<-spval
    for(i in (m-1):1){
      tmp[i]<-min(tmp[i+1],min((m/i)*spval[i],1,na.rm=TRUE),na.rm=TRUE)
      if(is.na(spval[i])) tmp[i]<-NA
    }
    adjp[,"BH"]<-tmp
  }

  if(is.element("BY",proc))
  {
    tmp<-spval
    a<-sum(1/(1:m))
    tmp[m]<-min(a*spval[m], 1) 
    for(i in (m-1):1){
      tmp[i]<-min(tmp[i+1],min((m*a/i)*spval[i],1,na.rm=TRUE),na.rm=TRUE)
      if(is.na(spval[i])) tmp[i]<-NA
    }
    adjp[,"BY"]<-tmp
  }

  if(is.element("ABH",proc))
  {
    # First obtain estimate of h_0, the number of true null hypotheses.
    tmp<-spval
    h0.m <- rep(0,m)
    for(k in 1:m){
      h0.m[k] <- (m+1-k)/(1-spval[k])
    }
    grab <- min(which(diff(h0.m,na.rm=TRUE)>0),na.rm=TRUE)
    h0.ABH <- ceiling(min(h0.m[grab],m))
    # Now apply BH procedure with adaptive correction.  
    for(i in (m-1):1){
      tmp[i]<-min(tmp[i+1],min((m/i)*spval[i],1,na.rm=TRUE),na.rm=TRUE)
      if(is.na(spval[i])) tmp[i]<-NA
    }
    adjp[,"ABH"]<-tmp*h0.ABH/m
  }

  list(adjp=adjp,index=index,h0.ABH=h0.ABH[1],h0.TSBH=h0.TSBH[1:length(alpha)])
}


###########################################################################

mt.reject<-function(adjp,alpha)
{
  which<-adjp<=alpha[1]
  dimnames(which)<-dimnames(adjp)

  if(is.matrix(adjp))
  {
    r<-matrix(0,length(alpha),ncol(adjp))
    for(i in 1:length(alpha))
        r[i,] <- colSums(adjp<=alpha[i])
    dimnames(r)<-list(alpha,dimnames(adjp)[[2]])
  }

  if(!is.matrix(adjp))
  {
    r<-rep(0,length(alpha))
    for(i in 1:length(alpha))
      r[i]<-sum(adjp<=alpha[i])
  }

  list(r=r,which=which)
}


###########################################################################

#need ... arg to legend to use with ... in mt.plot
mt.legend<-function(x, y = NULL, legend, fill = NULL, col = "black", lty, 
    lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"), 
    pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, xjust = 0, 
    yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), 
    text.width = NULL, text.col = par("col"), merge = do.lines && 
        has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE,...) 
{
    if (missing(legend) && !missing(y) && (is.character(y) || 
        is.expression(y))) {
        legend <- y
        y <- NULL
    }
    mfill <- !missing(fill) || !missing(density)
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) 
        stop("invalid coordinate lengths")
    xlog <- par("xlog")
    ylog <- par("ylog")
    rect2 <- function(left, top, dx, dy, density = NULL, angle, 
        ...) {
        r <- left + dx
        if (xlog) {
            left <- 10^left
            r <- 10^r
        }
        b <- top - dy
        if (ylog) {
            top <- 10^top
            b <- 10^b
        }
        rect(left, top, r, b, angle = angle, density = density, 
            ...)
    }
    segments2 <- function(x1, y1, dx, dy, ...) {
        x2 <- x1 + dx
        if (xlog) {
            x1 <- 10^x1
            x2 <- 10^x2
        }
        y2 <- y1 + dy
        if (ylog) {
            y1 <- 10^y1
            y2 <- 10^y2
        }
        segments(x1, y1, x2, y2, ...)
    }
    points2 <- function(x, y, ...) {
        if (xlog) 
            x <- 10^x
        if (ylog) 
            y <- 10^y
        points(x, y, ...)
    }
    text2 <- function(x, y, ...) {
        if (xlog) 
            x <- 10^x
        if (ylog) 
            y <- 10^y
        text(x, y, ...)
    }
    if (trace) 
        catn <- function(...) do.call(cat, c(lapply(list(...), 
            formatC), list("\n")))
    cin <- par("cin")
    Cex <- cex * par("cex")
    if (is.null(text.width)) 
        text.width <- max(strwidth(legend, units = "user", cex = cex))
    else if (!is.numeric(text.width) || text.width < 0) 
        stop("text.width must be numeric, >= 0")
    xc <- Cex * xinch(cin[1], warn.log = FALSE)
    yc <- Cex * yinch(cin[2], warn.log = FALSE)
    xchar <- xc
    yextra <- yc * (y.intersp - 1)
    ymax <- max(yc, strheight(legend, units = "user", cex = cex))
    ychar <- yextra + ymax
    if (trace) 
        catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, 
            ychar))
    if (mfill) {
        xbox <- xc * 0.8
        ybox <- yc * 0.5
        dx.fill <- xbox
    }
    do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
        0))) || !missing(lwd)
    n.leg <- if (is.call(legend)) 
        1
    else length(legend)
    n.legpercol <- if (horiz) {
        if (ncol != 1) 
            warning("horizontal specification overrides: Number of columns := ", 
                n.leg)
        ncol <- n.leg
        1
    }
    else ceiling(n.leg/ncol)
    if (has.pch <- !missing(pch) && length(pch) > 0) {
        if (is.character(pch) && !is.na(pch[1]) && nchar(pch[1]) > 
            1) {
            if (length(pch) > 1) 
                warning("Not using pch[2..] since pch[1] has multiple chars")
            np <- nchar(pch[1])
            pch <- substr(rep.int(pch[1], np), 1:np, 1:np)
        }
        if (!merge) 
            dx.pch <- x.intersp/2 * xchar
    }
    x.off <- if (merge) 
        -0.7
    else 0
    if (xlog) 
        x <- log10(x)
    if (ylog) 
        y <- log10(y)
    if (nx == 2) {
        x <- sort(x)
        y <- sort(y)
        left <- x[1]
        top <- y[2]
        w <- diff(x)
        h <- diff(y)
        w0 <- w/ncol
        x <- mean(x)
        y <- mean(y)
        if (missing(xjust)) 
            xjust <- 0.5
        if (missing(yjust)) 
            yjust <- 0.5
    }
    else {
        h <- n.legpercol * ychar + yc
        w0 <- text.width + (x.intersp + 1) * xchar
        if (mfill) 
            w0 <- w0 + dx.fill
        if (has.pch && !merge) 
            w0 <- w0 + dx.pch
        if (do.lines) 
            w0 <- w0 + (2 + x.off) * xchar
        w <- ncol * w0 + 0.5 * xchar
        left <- x - xjust * w
        top <- y + (1 - yjust) * h
    }
    if (plot && bty != "n") {
        if (trace) 
            catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
                h, ", ...)", sep = "")
        rect2(left, top, dx = w, dy = h, col = bg, density = NULL)
    }
    xt <- left + xchar + (w0 * rep.int(0:(ncol - 1), rep.int(n.legpercol, 
        ncol)))[1:n.leg]
    yt <- top - (rep.int(1:n.legpercol, ncol)[1:n.leg] - 1) * 
        ychar - 0.5 * yextra - ymax
    if (mfill) {
        if (plot) {
            fill <- rep(fill, length.out = n.leg)
            rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
                col = fill, density = density, angle = angle, 
                border = "black")
        }
        xt <- xt + dx.fill
    }
    if (plot && (has.pch || do.lines)) 
        col <- rep(col, length.out = n.leg)
    if (missing(lwd)) 
        lwd <- par("lwd")
    if (do.lines) {
        seg.len <- 2
        if (missing(lty)) 
            lty <- 1
        ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
        lty <- rep(lty, length.out = n.leg)
        lwd <- rep(lwd, length.out = n.leg)
        if (trace) 
            catn("  segments2(", xt[ok.l] + x.off * xchar, ",", 
                yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
        if (plot) 
            segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * 
                xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l], 
                col = col[ok.l])
        xt <- xt + (seg.len + x.off) * xchar
    }
    if (has.pch) {
        pch <- rep(pch, length.out = n.leg)
        pt.bg <- rep(pt.bg, length.out = n.leg)
        pt.cex <- rep(pt.cex, length.out = n.leg)
        pt.lwd <- rep(pt.lwd, length.out = n.leg)
        ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
        x1 <- (if (merge) 
            xt - (seg.len/2) * xchar
        else xt)[ok]
        y1 <- yt[ok]
        if (trace) 
            catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
                ", ...)")
        if (plot) 
            points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
                bg = pt.bg[ok], lwd = pt.lwd[ok])
        if (!merge) 
            xt <- xt + dx.pch
    }
    xt <- xt + x.intersp * xchar
    if (plot) 
        text2(xt, yt, labels = legend, adj = adj, cex = cex, 
            col = text.col)
    invisible(list(rect = list(w = w, h = h, left = left, top = top), 
        text = list(x = xt, y = yt)))
}

mt.plot<-function(adjp,teststat, plottype="rvsa",logscale=FALSE,
                  alpha=seq(0,1,length=100), proc="",leg=c(0,0),...)
{
  m<-nrow(adjp)
  n<-ncol(adjp)
  a<-length(alpha)

  if(plottype=="rvsa")
  {
    r<-mt.reject(adjp,alpha)$r
    matplot(alpha,r,xlab="Type I error rate",
            ylab="Number of rejected hypotheses", type="l", ...)
    mt.legend(leg[1],leg[2],proc,...)
  }

  if(plottype=="pvsr")
  {
    spval<-apply(adjp,2,sort)
    matplot(1:m,spval,xlab="Number of rejected hypotheses",
            ylab="Sorted adjusted p-values", type="l", ...)
    mt.legend(leg[1],leg[2],proc,...)
  }

  if(plottype=="pvst")
  {
    if(!logscale)
      matplot(teststat,adjp,xlab="Test statistics",
              ylab="Adjusted p-values", type="p", ...)
    if(logscale)
      matplot(teststat,-log(adjp,10),xlab="Test statistics",
              ylab="-log(adjusted p-values,10)", type="p", ...)
    mt.legend(leg[1],leg[2],proc,...)
  }
  if(plottype=="pvsi")
  {
    if(!logscale)
      matplot(1:m,adjp,xlab="index",ylab="Adjusted p-values", type="l", ...)
    if(logscale)
      matplot(1:m,-log(adjp,10),xlab="index",
              ylab="-log(adjusted p-values,10)", type="l", ...)
    mt.legend(leg[1],leg[2],proc,...)
  }
}



