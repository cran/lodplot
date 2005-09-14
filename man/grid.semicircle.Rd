\name{grid.semicircle}
\alias{grid.semicircle}
\title{ Draws a semicircle }
\description{
Using grid graphical commands, draws a semicircle (including the base).
}
\usage{
grid.semicircle(base.x, base.y, base.length, height = base.length, 
  side = 1, orientation = NULL, col = NULL)
}
\arguments{
  \item{base.x}{is X coordinate of base }
  \item{base.y}{is Y coordinate of base }
  \item{base.length}{is length of base of semicircle}
  \item{height}{is height of semicircle}
  \item{side}{is orientation of semicircle, one of four cardinal directions }
  \item{orientation}{is orientation, in degrees}
  \item{col}{is line colour used to draw the semicircle}
}
\details{
Draw a semicircle either oriented to the four cardinal points (\code{side}
1=below 2=left 3=above, 4=right), in which case \code{base.x} and
\code{base.y} determine the bottom left point of base, or an angle,
when \code{base.x} and \code{base.y} are the midpoint of the base.
}
\value{
A grid graphical object.
}
\author{ David L Duffy }
\keyword{ iplot }
