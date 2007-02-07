\name{chromosome.viewsequence}
\alias{chromosome.viewsequence}
\title{  Produce a viewport with the chromosome along the top
and a pendant lod score plot
}
\description{
Pretty plot of location score versus sequence map position on one 
(human) chromosome
}
\usage{
chromosome.viewsequence(x, chrom, statistic = "lod", pheno.names = NULL, 
  min.stat = 0, max.stat = 4, col = 1:6, lwd = 2, lty = 1, 
  hpos = 0.85, width = 0.05, chromname.cex=1.5, 
  units = "bp", bands = "major", xticdist=5e7, 
  show.y.axis = FALSE, new = FALSE, ...)
}
\arguments{
  \item{x}{is a data.frame containing variables \code{chr} (indicating the
chromosome 1..X), \code{pos} (the map position, in base pairs), and the
statistics to be plotted.}
  \item{chrom}{is the chromosome to be plotted.}
  \item{statistic}{is the vector of test statistics to be plotted, 
defaulting to "lod".}
  \item{pheno.names}{gives a long name for the statistics, defaulting to 
the statistic name.}
  \item{min.stat}{is the minimum plottable value of the statistic.}
  \item{max.stat}{is the maximum plottable value of the statistic.}
  \item{col}{is a vector of colours for the lod curve.}
  \item{lty}{is a vector of line types for the lod curve.}
  \item{lwd}{is a vector of line widths for the lod curve.}
  \item{chromname.cex}{is relative font size for the chromosome label.}
  \item{units}{are the map units, defaulting to "bp".}
  \item{hpos}{is the Y coordinate for the chromosome ideogram }
  \item{width}{is the width of the chromosome ideogram }
  \item{bands}{is which chromosomal bands to display on the ideogram.}
  \item{xticdist}{is the spacing of the X axis tickmarks, defaulting to 5e7.}
  \item{show.y.axis}{determines whether the Y axis should be shown.}
  \item{new}{indicates whether addition to existing plot.}
  \item{\dots}{are other graphical parameters to be passed to plot.default.}
}
\details{
Given a set of genetic location scores and chromosomal positions,
\code{chromosome.viewsequence} produces a plot for chromosome \code{chrom}.
A chromosome ideogram is drawn at the top of the plot, with band locations
expressed on the same scale as the map positions (bp).
}
\value{
A grid graphical object.
}
\author{ David L Duffy }
\keyword{ hplot }
