"paint.chromosome" <-
function(chrom, pos=0, units="cM", width=0.4, bands="major") {
#
# base graphics semicircle 
  semicircle <- function(base.x, base.y, base.length,
                         height=base.length, side=1, orientation=NULL, col=NULL) {
    radius<-base.length/2
    x<-radius*seq(-1,1,length=40)
    y<-height/radius*sqrt(radius^2-x^2)
    if (is.null(orientation)) {
      co<-as.integer(cos(pi*(3-side)/2))
      so<-as.integer(sin(pi*(3-side)/2))
    }else{
      co<-cos(orientation)
      so<-sin(orientation)
    }
    tx<-co*x - so*y 
    ty<-so*x + co*y
    if (is.null(orientation)) {
      if (side==1 || side==3) {
        base.x<-base.x+radius
      }else if (side==2 || side==4) {
        base.y<-base.y+radius
      }
    }
    x<-base.x+tx
    y<-base.y+ty
    polygon(x,y,col=col)
  }
  data(chrom.bands)
  chromdata<-subset(chrom.bands, chrom.bands$chr==chrom)
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
  rect(chromdata$cM.top[idx],pos,chromdata$cM.bot[idx],pos-width, col=bandcol[idx])
  semicircle(chromdata$cM.bot[1], pos-width, width,
             chromdata$cM.bot[1]-chromdata$cM.top[1], 2, col=bandcol[1])
  semicircle(chromdata$cM.top[n], pos-width, width, 
             chromdata$cM.bot[n]-chromdata$cM.top[n], 4, col=bandcol[n])
  semicircle(chromdata$cM.top[centromere], pos-width, width,
             chromdata$cM.bot[centromere]-chromdata$cM.top[centromere], 
             4, col=bandcol[centromere])
  semicircle(chromdata$cM.bot[centromere+1], pos-width, width, 
             chromdata$cM.bot[centromere+1]-chromdata$cM.top[centromere+1], 
             2, col=bandcol[centromere+1])
  points(chromdata$cM.top[centromere], pos-0.5*width, col="black", cex=3, pch=16)
  points(chromdata$cM.top[centromere], pos-0.5*width, col="white", cex=3, pch=20)
}

