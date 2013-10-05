\name{BG}
\alias{BG}
\title{
  BG test data set
}
\description{
  The BG data set has 10 descriptors evaluated by comparison of a
  trial product with a reference product on an ordinal scale (-2, -1,
  0, 1, 2). Thirty tasters made two comparisons (ref-ref,
  ref-trial) on all products yielding 60 observations in  total.
}
\usage{
data(BG)
}
\format{
\describe{
  \item{\code{Taster}}{
    factor with 11 levels.
  }
  \item{\code{Samples}}{
    factor with 2 levels: \code{"565", "987"} indicating the
    two different comparisons made: ref-ref, ref-trial.
  }
  \item{\code{Yellow color, Shine, Granules, Thickness, Stringy}}{
    Descriptor variables; intergers taking values between -2 and 2.
  }
  \item{\code{Smooth, Sticky in mouth, Acidity, Global Dairy, Off notes}}{
    Descriptor variables; intergers taking values between -2 and 2.
  }
}}
\seealso{
  \code{\link{opair}}.
}
\source{
  Data are produced by Danone Research.
}
%% \references{
%% }
\examples{

## get vector of descriptor/attribute names:
(attr.names <- names(BG)[3:ncol(BG)])

## Fit Thurstonian ordinal paired comparison model:
(fit <- opair(BG[attr.names], BG$Samples, d.equiv=1,
conf.level=.90))

## plot d-primes:
plot(fit)
plot(fit, type=2)

}
\keyword{datasets}
