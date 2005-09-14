\name{halfsibscan}
\alias{halfsibscan}
\title{ Simulation of affected half-sib pair genome scan
}
\description{
Generate results from an affected half-sib pair genome scan
}
\usage{
halfsibscan(N = 100, grid = 1)
}
\arguments{
  \item{N}{is the number of affected half-sib pair families}
  \item{grid}{the distance between the (perfectly informative) scan markers}
}
\details{
Carried out the appropriate random walk on each chromosome to
generate lod scores that would result if there were no
trait loci segregating in the pedigrees sampled.
}
\value{
A data frame, containing: 
  \item{\code{chr}}{the chromosome (1..22, X)}
  \item{\code{pos}}{the genetic map position on the chromosome}
  \item{\code{lod}}{the simulated lod score.}
}
\author{ David L Duffy }
\seealso{ plot.scan }
\examples{
#
# Plot one simulated scan
#
  plot.scan(halfsibscan(), with.X=TRUE)
}
\keyword{ datagen }
