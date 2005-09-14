#
# random walk in affected half-sibs
#
halfsibscan <- function(N=100, grid=1) {
  chr <- c(1:22, "X")
  chr.length <- c(285, 265, 225, 210, 220, 190, 190, 175, 165, 175,
                  160, 175, 130, 120, 130, 130, 135, 135, 110, 100,
                  70, 80, 185)
  halfsibZ <- function(N, grid, chr, length) {
    phi <- 0.5*(1-exp(-0.04*grid))
    npoints <- round(length/grid)+1
    points <- seq(from=0, to=length, length.out=npoints)
    y <- double(length=npoints)
    y[1] <- sum(sample(0:1, replace=TRUE, size=N))
    for(i in 2:length) {
      y[i] <- y[i-1] - rbinom(1, y[i-1], phi) + rbinom(1, N-y[i-1], phi)
    }
    z <- (2*y-N)/sqrt(N)
    z[z<0] <- 0
    data.frame(chr=rep(chr,npoints), pos=points, lod=z^2/(2*log(10)))
  }
  res <- NULL
  for(i in 1:23) {
    res <- rbind(res, halfsibZ(N=N, grid=grid, chr=chr[i], length=chr.length[i]))
  }
  res
}
