\name{save.opair}
\alias{save.opair}
\title{
  Save d-prime from Thurstonian model fit to ordinal paired comparisons
}
\description{
  This function saves the coefficient table (including d-prime values)
  to a csv file from an \code{opair} fit.
}
\usage{
  save.opair(fit, file=stop("'file' must be specified"), sep=",",
             dec=".", ...)
}
\arguments{
  \item{fit}{
    a \code{opair} object produced by the \code{opair} function.
  }
  \item{file}{
    a character vector naming the resulting file. The \code{".csv"}
    ending is automatically added, thus using, say, \code{file="name"}
    produces the file \code{"name.csv"}.
  }
  \item{sep}{
    the column separator, cf. \code{\link{write.table}}
  }
  \item{dec}{
    the decimal separator, cf. \code{\link{write.table}}
  }
  \item{\dots}{
    additional arguments passed on to \code{\link{write.table}}.
  }
}
\value{
  No return value - the function is used for its side effects (the
  \code{"csv"} file).
}
\details{
  The function produces a \code{"csv"} file using

  \code{write.table(coef(fit), file=paste0("name", ".csv"), quote=FALSE,
    sep=sep, row.names=FALSE, dec=dec, ...}.
}
%% \references{
%% }
\seealso{
  \code{\link{opair}}, \code{\link{write.table}}
}
\examples{

## Fit opair model to the NV data set:
(attr.names <- names(NV)[4:ncol(NV)])
fit <- opair(NV[, attr.names], NV$Samples,
             abbreviate.names=FALSE)

## Save the coefficient table from a opair fit to a csv file:
\dontrun{
save.opair(fit, "myDprimes", sep=";", dec=",")
}
## The resulting file is called "myDprimes.csv" and it is placed in the
## current working directory:
getwd()

}
\author{Rune Haubo B Christensen}
\keyword{file}
