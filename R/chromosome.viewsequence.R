"chromosome.viewsequence" <-
function(x, chrom, statistic="lod", 
                            pheno.names=NULL, 
                            min.stat=0, max.stat=4, 
                            col=1:6, lwd=2, lty=1, 
                            hpos=0.85, width=0.05, 
                            chromname.cex=1.5,   
                            units="bp", bands="major", xticdist=5e7,
                            show.y.axis=FALSE, new=FALSE, ...) {
# chromosome
  require(grid)
  data(chrom.bands)
  if (new) grid.newpage()
  chromdata<-subset(chrom.bands, chrom.bands$chr==chrom)
  if (nrow(chromdata)==0) stop(paste("No band data for chromosome ",chrom,"!",sep=""))
  lc<-nchar(chromdata$band)
  sel<-!(substr(chromdata$band,lc,lc) %in% letters)
  if (bands!="major") sel<-!sel
  chromdata<-chromdata[sel,]
  rm(lc,sel)
  bandcol<-gray(c(0.4,0.6,0.8,0.8,0.85))[match(chromdata$stain, 
                                          c("acen","gneg", "gpos", "gvar", "stalk"))]
  n<-nrow(chromdata)
  centromere<-which(chromdata$arm[-n]!=chromdata$arm[-1])
  idx<-c(2:(centromere-1), (centromere+2):(n-1))
  pushViewport(viewport(xscale=c(chromdata$bases.top[1]-5,chromdata$bases.bot[n]+5), 
                        yscale=c(0,1),
                        clip="on"))
  grid.rect(x=chromdata$bases.top[idx],y=hpos,
            width=chromdata$bases.bot[idx]-chromdata$bases.top[idx],
            height=width,
            just=c("left","top"),
            default.units="native", gp=gpar(fill=bandcol[idx]))
  grid.semicircle(chromdata$bases.bot[1], hpos-width, width,
                  chromdata$bases.bot[1]-chromdata$bases.top[1], 2, col=bandcol[1])
  grid.semicircle(chromdata$bases.top[n], hpos-width, width, 
                  chromdata$bases.bot[n]-chromdata$bases.top[n], 4, col=bandcol[n])
  grid.semicircle(chromdata$bases.top[centromere], hpos-width, width,
                  chromdata$bases.bot[centromere]-chromdata$bases.top[centromere], 
                  4, col=bandcol[centromere])
  grid.semicircle(chromdata$bases.bot[centromere+1], hpos-width, width, 
                  chromdata$bases.bot[centromere+1]-chromdata$bases.top[centromere+1], 
                  2, col=bandcol[centromere+1])
  grid.points(unit(chromdata$bases.bot[centromere],"native"), 
              unit(hpos-0.5*width,"native"),
              size=unit(1.5,"char"), pch=20, gp=gpar(col="white"))
  grid.points(unit(chromdata$bases.bot[centromere],"native"), 
              unit(hpos-0.5*width,"native"),
              size=unit(0.5,"char"), pch=20, gp=gpar(col="black"))
  grid.text(chrom,
            unit(0.5,"npc"),
            unit(hpos+2*width,"native"), gp=gpar(cex=chromname.cex))
# stat curve
  pos<-x$pos[x$chr %in% chrom]
  stat<-x[x$chr %in% chrom, statistic]
  hi.stat<-max(unlist(c(max.stat,stat)),na.rm=TRUE)
  if (sum(x$chr %in% chrom)>0) {
    pushViewport(viewport(x=unit(0,"native") ,y=unit(1,"lines"),
                          width=unit(max(pos, na.rm=TRUE),"native"),
                          height=unit(0.95, "npc"),
                          just=c("left","bottom"), 
                          xscale=c(0,max(pos,na.rm=TRUE)), 
                          yscale=c(min.stat,hi.stat),
                          clip="off"))
    grid.rect()
    if (is.vector(stat)) {
      grid.lines(pos, stat, default.units="native", 
      gp=gpar(col=col, lwd=lwd, lty=lty, ...))
    }else{
      if (length(col)<ncol(stat)) {
        col <- rep(col, length.out = ncol(stat))
      }
      if (length(lty)<ncol(stat)) {
        lty <- rep(lty, length.out = ncol(stat))
      }
      for(j in 1:ncol(stat)) {
        grid.lines(pos,stat[,j], default.units="native", 
        gp=gpar(col=col[j], lwd=lwd, lty=lty[j], ...))
      }
    }
    yticks <- seq(min.stat, max.stat,1)
    my.xtics(at=seq(0,max(pos,na.rm=TRUE),xticdist),length=0.25)
    grid.grill(h=yticks, v=0, 
               default.units="native",gp=gpar(lty=3))
    if (show.y.axis || new) {
      grid.text(yticks,
              unit(rep(0.0,4),"npc")-unit(rep(0.5,4),"lines"),
              unit(0:3,"native"), gp=gpar(cex=0.75))
    }
    popViewport()
  }
  popViewport()
}

