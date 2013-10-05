\name{plot}
\alias{plot.opair}
\title{
  Plot d-prime from Thurstonian model fit to ordinal paired comparisons
}
\description{
  This function provides a dot-plot or a bar-plot of d-prime values from
  the fit of a Thurstonian model to ordinal paired comparisons.
}
\usage{
  \method{plot}{opair}(x, type=1, ...)
}
\arguments{
  \item{x}{
    an \code{opair} object produced by the \code{opair} function.
  }
  \item{type}{
    a numerical scalar: \code{1} produces a dot-plot (using
    \code{dotchart}) and \code{2} produces a bar-plot (using
    \code{plot} with \code{type="h"}).
  }
  \item{\dots}{
    additional arguments passed on to \code{dotchart} or \code{plot}.
  }
}
\value{
  No return value - the function is used for its side effects (the
  plots).
}
%% \details{
%% }
%% \references{
%% }
\seealso{
  \code{\link{opair}}
}
\examples{

## Fit opair model to the NV data set:
(attr.names <- names(NV)[4:ncol(NV)])
(fit <- opair(NV[, attr.names], NV$Samples))

## dot plot:
plot(fit) ## or plot(fit, type=1)

## bar plot:
plot(fit, type=2)

}
\author{Rune Haubo B Christensen}
\keyword{hplot}
