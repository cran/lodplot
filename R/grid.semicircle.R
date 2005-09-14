"grid.semicircle" <-
function(base.x, base.y, base.length,
                            height=base.length, side=1, 
                            orientation=NULL, col=NULL) {
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
  grid.polygon(x,y,default.units="native", gp=gpar(fill=col))
}

