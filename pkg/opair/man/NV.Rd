\name{NV}
\alias{NV}
\title{
  NV test data set
}
\description{
  The NV data set has 12 descriptors evaluated by comparison of two
  trial products with one reference product on an ordinal scale (-2, -1,
  0, 1, 2). Eleven tasters made all three comparisons (ref-ref,
  ref-trial1, ref-trial2) on all products yielding 33 observations in
  total.
}
\usage{
data(NV)
}
\format{
\describe{
  \item{\code{Taster}}{
    factor with 11 levels.
  }
  \item{\code{Samples}}{
    factor with 3 levels: \code{"432", "568", "841"} indicating the
    three different comparisons made: ref-ref, ref-trial1, ref-trial2.
  }
  \item{\code{Red color, Thickness, Fruit quantity, Smooth,
      Compact, Sweetness}}{
    Descriptor variables; intergers taking values between -2 and 2.
  }
  \item{\code{Acidity, Fruit, Global Dairy, Off notes,
      Fruit intensity, Lasting sweet}}{
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
(attr.names <- names(NV)[3:ncol(NV)])

## Fit Thurstonian ordinal paired comparison model:
(fit <- opair(NV[attr.names], NV$Samples, d.equiv=1))

## plot d-primes:
plot(fit)
plot(fit, type=2)

}
\keyword{datasets}
