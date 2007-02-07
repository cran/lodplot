"plot.scan" <-
function(x, type="layout", statistic="lod", 
                      with.X=TRUE, min.stat=0, max.stat=4, pheno.names=NULL,
                      units="cM", col=1:6, lty=1, lwd=2,
                      chromname.cex=0.9, ...) {
  require(grid)
  if (with.X) {
    n.chrom<-23
    nchrom <- c(3,4,7,9)
    chroms<-c(as.character(1:22),"X")
    chrom.rel.size<- list(r1=c(5, 5, 5), 
                          r2=c(4, 4, 4, 3), 
                          r3=c(3, 2, 2, 2, 2, 2, 2),
                          r4=c(2, 2, 2, 2, 1, 1, 1, 1, 3))
  }else{
    n.chrom<-22
    nchrom <- c(3,4,7,8)
    chroms<-as.character(1:22)
    chrom.rel.size<- list(r1=c(5, 5, 5), 
                          r2=c(4, 4, 4, 3), 
                          r3=c(3, 2, 2, 2, 2, 2, 2),
                          r4=c(2, 2, 2, 2, 2, 2, 2, 2))
  }
  type<-match.arg(type,c("layout","overwrite","linear","histogram"))
  if (is.null(pheno.names)) {
    pheno.names <- statistic
  }
  if (type=="layout") {
    grid.newpage()
    chr.plot <- list(
      r1=viewport(1/2, 7/8, width=0.9, height=0.24, name="r1", 
                layout=grid.layout(nr=1, nc=3, 
                widths=unit(chrom.rel.size$r1, "null"))), 
      r2=viewport(1/2, 5/8, width=0.9, height=0.24, name="r2", 
                layout=grid.layout(nr=1, nc=4, 
                widths=unit(chrom.rel.size$r2, "null"))),
      r3=viewport(1/2, 3/8, width=0.9, height=0.24, name="r3", 
                layout=grid.layout(nr=1, nc=7, 
                widths=unit(chrom.rel.size$r3, "null"))), 
      r4=viewport(1/2, 1/8, width=0.9, height=0.24, name="r4", 
                   layout=grid.layout(nr=1, nc=length(chrom.rel.size$r4), 
                            widths=unit(chrom.rel.size$r4, "null"))))
    i.chr <- 0
    for(g in seq(along=chr.plot)) {
      pushViewport(chr.plot[[g]])
      if (units=="bp") {
        for(i in 1:nchrom[g]) {
          i.chr <- i.chr + 1
          pushViewport(viewport(layout.pos.col=i, layout.pos.row=1))
          chromosome.viewsequence(x,chroms[i.chr], statistic,col=col, 
                                  min.stat=min.stat, max.stat=max.stat, 
                                  lwd=lwd, lty=lty, 
                                  chromname.cex=chromname.cex, 
                                  show.y.axis=(i==1))
          popViewport()
        }
      }else{
        for(i in 1:nchrom[g]) {
          i.chr <- i.chr + 1
          pushViewport(viewport(layout.pos.col=i, layout.pos.row=1))
          chromosome.viewlinkage(x,chroms[i.chr],statistic,col=col, 
                                 min.stat=min.stat, max.stat=max.stat, 
                                 lwd=lwd, lty=lty, 
                                 chromname.cex=chromname.cex, 
                                 show.y.axis=(i==1))
          popViewport()
        }
      }
      popViewport()
    }
  }else if (type=="linear") {
    n.chrom<-length(unique(x$chr))
    statlod<-x[, statistic]
    hi.lod<-max(unlist(c(max.stat,lod)),na.rm=TRUE)
    pos<-x$pos
    d<-diff(pos)
    d[d<0]<-0
    pos<-c(0,cumsum(d))
    tic.pos<-0
    tic.lab<-"0"
    tic.lab.last<-0
    cur.chr<-x$chr[1]
    for (j in 1:length(pos)) {
      if (x$chr[j]!=cur.chr) {
        tic.pos<-c(tic.pos, pos[j])
        tic.lab<-c(tic.lab,"")
        tic.lab.last<-0
        cur.chr<-x$chr[j]
      }else if (pos[j]>(tic.pos[length(tic.pos)]+50)) {
        tic.pos<-c(tic.pos, tic.pos[length(tic.pos)]+50)
        tic.lab.last<-tic.lab.last+50
        tic.lab<-c(tic.lab,as.character(tic.lab.last))
      }
    }
    mid.pos<-tapply(pos, x$chr, mean, na.rm=TRUE)
    if (is.null(ncol(lod))) {
      npheno<-1
      plot(pos, lod, t="l",
           ylim=c(min.stat,hi.lod), axes=FALSE,
           xlab=paste("Genome scan position (",units,")",sep=""), 
           ylab="lod score", ...)
    }else{
      npheno <- ncol(lod)
      matplot(pos, lod, t="l",
           ylim=c(min.stat,hi.lod), axes=FALSE,
           col=col, lwd=lwd, 
           xlab=paste("Genome scan position (",units,")",sep=""), 
           ylab="lod score", ...)
    }
    box()
    if (n.chrom>15) {
      axis(1)
    }else{
      axis(1, at=tic.pos, labels=tic.lab, cex.axis=min(1,50/length(tic.pos)))
    }
    axis(2)
    abline(v=pos[d==0], lty=3)
    text(mid.pos, 3.5, names(mid.pos), cex=chromname.cex)
    abline(h=3, lty=3)
    abline(h=0, lty=3)
    if (npheno>1) {
      legend(0, max.stat+0.1, lwd=rep(2,npheno), 
             bty="o",  bg="white", 
             col=c("black","red", "blue"),
             legend=pheno.names, horiz=TRUE)
    }
  }else if (type=="histogram") {
    hist(x[,statistic],20,xlab="lod score",
         main="Distribution of lod scores across scan", ...)
  }else{
    lod<-x[, statistic[1]]
    hi.lod<-max(unlist(c(max.stat,lod)),na.rm=TRUE)
    idx<-as.numeric(as.character(x$chr))==1
    plot(x$pos[idx], lod[idx], t="l",
         ylim=c(min.stat,hi.lod), 
         xlab=paste("Map position (i",units,")",sep=""), ylab="lod score", ...)
    for (i in unique(x$chr)) {
      lty<-1
      idx<-(x$chr==i)
      if (sum(idx)>0) {
        lty<-as.numeric(as.character(i))
        lines(x$pos[idx], lod[idx], lty=lty)
      }
    }
  }
}

