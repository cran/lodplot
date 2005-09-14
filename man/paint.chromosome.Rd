\name{paint.chromosome}
\alias{paint.chromosome}
\title{ Draw a chromosome ideogram }
\description{
Produce an ideogram of a specified  human chromosome, using base graphics.
Usually for addition to a plot of linkage analysis results, so the
chromosome length is in genetic map units (usually cM).
}
\usage{
paint.chromosome(chrom, pos = 0, units = "cM", width = 0.4, bands = "major")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{is the chromosome (1..X).}
  \item{pos}{is the Y position.}
  \item{units}{is the user X axis units, usually "cM" or "M".}
  \item{width}{is the width of the ideogram.}
  \item{bands}{is which chromosomal bands to display on the ideogram.}
}
\value{
A base graphics plot.
}
\author{ David L Duffy }
\keyword{ aplot }
