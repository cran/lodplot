"chromosome.view" <-
function(x, chrom, statistic="lod", 
                            pheno.names=NULL, 
                            min.lod=0, max.lod=4, 
                            col=1:6, lwd=2, lty=1, 
                            hpos=0.85, width=0.05, 
                            chromname.cex=1.5,   
                            units="cM", bands="major", 
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
  pushViewport(viewport(xscale=c(chromdata$cM.top[1]-5,chromdata$cM.bot[n]+5), 
                        yscale=c(0,1),
                        clip="on"))
  grid.rect(x=chromdata$cM.top[idx],y=hpos,
            width=chromdata$cM.bot[idx]-chromdata$cM.top[idx],
            height=width,
            just=c("left","top"),
            default.units="native", gp=gpar(fill=bandcol[idx]))
  grid.semicircle(chromdata$cM.bot[1], hpos-width, width,
                  chromdata$cM.bot[1]-chromdata$cM.top[1], 2, col=bandcol[1])
  grid.semicircle(chromdata$cM.top[n], hpos-width, width, 
                  chromdata$cM.bot[n]-chromdata$cM.top[n], 4, col=bandcol[n])
  grid.semicircle(chromdata$cM.top[centromere], hpos-width, width,
                  chromdata$cM.bot[centromere]-chromdata$cM.top[centromere], 
                  4, col=bandcol[centromere])
  grid.semicircle(chromdata$cM.bot[centromere+1], hpos-width, width, 
                  chromdata$cM.bot[centromere+1]-chromdata$cM.top[centromere+1], 
                  2, col=bandcol[centromere+1])
  grid.points(unit(chromdata$cM.bot[centromere],"native"), 
              unit(hpos-0.5*width,"native"),
              size=unit(1.5,"char"), pch=20, gp=gpar(col="white"))
  grid.points(unit(chromdata$cM.bot[centromere],"native"), 
              unit(hpos-0.5*width,"native"),
              size=unit(0.5,"char"), pch=20, gp=gpar(col="black"))
  grid.text(chrom,
            unit(0.5,"npc"),
            unit(hpos+2*width,"native"), gp=gpar(cex=chromname.cex))
# lod curve
  pos<-x$pos[x$chr %in% chrom]
  lod<-x[x$chr %in% chrom, statistic]
  hi.lod<-max(unlist(c(max.lod,lod)),na.rm=TRUE)
  if (sum(x$chr %in% chrom)>0) {
    pushViewport(viewport(x=unit(0,"native") ,y=unit(1,"lines"),
                          width=unit(max(pos, na.rm=TRUE),"native"),
                          height=unit(0.95, "npc"),
                          just=c("left","bottom"), 
                          xscale=c(0,max(pos,na.rm=TRUE)), 
                          yscale=c(min.lod,hi.lod),
                          clip="off"))
    grid.rect()
    if (is.vector(lod)) {
      grid.lines(pos, lod, default.units="native", 
      gp=gpar(col=col, lwd=lwd, lty=lty, ...))
    }else{
      if (length(col)<ncol(lod)) {
        col <- rep(col, length.out = ncol(lod))
      }
      if (length(lty)<ncol(lod)) {
        lty <- rep(lty, length.out = ncol(lod))
      }
      for(j in 1:ncol(lod)) {
        grid.lines(pos,lod[,j], default.units="native", 
        gp=gpar(col=col[j], lwd=lwd, lty=lty[j], ...))
      }
    }
    my.xtics(at=seq(0,max(pos,na.rm=TRUE),50),length=0.25)
    grid.grill(h=0:3, v=0, default.units="native",gp=gpar(lty=3))
    if (show.y.axis || new) {
      grid.text(0:3,
              unit(rep(0.0,4),"npc")-unit(rep(0.5,4),"lines"),
              unit(0:3,"native"), gp=gpar(cex=0.75))
    }
    popViewport()
  }
  popViewport()
}

