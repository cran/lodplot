\name{plot.scan}
\alias{plot.scan}
\title{ Plot results of a linkage genome scan }
\description{
Produces several types of plot of the lod score (or other test statistic)
versus genetic map position (cM) on each of 22 human autosomal chromosomes
along with the X chromosome if requested.  Decorates with a chromosome 
ideogram if appropriate.
}
\usage{
plot.scan(x, type = "layout", statistic = "lod", with.X = TRUE, 
  min.lod = 0, max.lod = 4, pheno.names = NULL, units = "cM", 
  col = 1:6, lty = 1, lwd = 2, chromname.cex=0.9, ...)
}
\arguments{
  \item{x}{is a data.frame containing variables \code{chr} (indicating the
chromosome 1..X), \code{pos} (the map position, usually in cM), and the
statistics to be plotted.}
  \item{type}{is the type of plot: "layout", "linear", "overwrite", or
  "histogram", defaulting to "layout".}
  \item{statistic}{is the test statistic to be plotted, defaulting to "lod".}
  \item{with.X}{flags whether the X chromosome is to be included in the plot.}
  \item{min.lod}{is the minimum plottable value of the statistic.}
  \item{max.lod}{is the maximum plottable value of the statistic.}
  \item{pheno.names}{gives a long name for the statistics, defaulting to 
the statistic name.}
  \item{units}{are the map units, defaulting to "cM".}
  \item{col}{is a vector of colours for the lod curve.}
  \item{lty}{is a vector of line types for the lod curve.}
  \item{lwd}{is a vector of line widths for the lod curve.}
  \item{chromname.cex}{is the relative font size for the chromosome label.}
  \item{\dots}{are other graphical parameters to be passed to plot.default.}
}
\details{
Given a set of genetic location scores and chromosomal positions,
\code{plot.scan} defaults to producing an array of plots ("layout"), one
per chromosome, with each scaled to the chromosome length.  A chromosome
ideogram is drawn at the top of each plot, with band locations expressed
on the same scale as the map positions (usually cM).

Alternatively, the "linear" option plots the entire genome on one graph, with
divisions between the chromosomes.  Chromosome ideograms are not added.
The "overwrite" option overwrites the curves for each chromosome onto one
graph, and the "histogram" option gives the distribution of lod scores over
the entire sample.

The function uses the grid library.  The ideogram band genetic map
positions were interpolated from the physical positions of the bands
estimated from the Build 35.1 human sequence for the NCBI Map Viewer.
}
\value{
A grid type graphical object.
}
\references{ 
Human chromosome ideograms.
ftp://ftp.ncbi.nlm.nih.gov/genomes/H\_sapiens/maps/mapview/BUILD.35.1/ideogram.gz.
Accessed September 2005.
}
\author{ David L Duffy }

\examples{
# Plot one simulated scan
  plot.scan(halfsibscan(), col="red", with.X=TRUE)
# plot a genome scan
\dontrun{x <-read.table("nfr.lod", h=T, na=".")}
\dontrun{plot.scan(x, statistic=c("flat.vc","flat.sqtl"), 
  col=c("red","blue"),with.X=TRUE)}
}
\keyword{ hplot }
