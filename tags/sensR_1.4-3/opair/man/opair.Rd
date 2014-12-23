\name{opair}
\alias{opair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Fit the Thurstonian model to ordinal paired comparisons
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  This function fits the Thurstonian model to ordinal paired comparisons
  for a series of descriptors on one or more products with maximum
  likelihood.
}
\usage{
opair(descriptors, products, d.equiv = 0.5, conf.level = 0.95,
abbreviate.names = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{descriptors}{
    %%     ~~Describe \code{descriptors} here~~
    a \code{data.frame} with the descriptor variables. Each variable
    is coerced into an ordered factor with levels and ordering given by
    the unique values of the variable. If a variable is already a factor
    or an ordered factor, it is first coerced into a numeric variable,
    so the resulting levels and possibly order might be different from
    the original ones.
  }
  \item{products}{
    %%     ~~Describe \code{products} here~~
    a factor-variable with two or more levels of the same length as
    descriptors. The first level is interpreted as the reference
    product.
  }
  \item{d.equiv}{
    %%     ~~Describe \code{d.equiv} here~~
    the point of equivalence for d-prime. The equivalence region is
    considered between \code{-d.equiv} and \code{d.equiv}.
  }
  \item{conf.level}{
    %%     ~~Describe \code{conf.level} here~~
    the confidence level for the confidence intervals for the d-prime
    estimates.
  }
  \item{abbreviate.names}{
    logical. Should descriptor and product names be abbreviated to
    preserve space in the printed output? Defaults to \code{FALSE}.
    }
  \item{\dots}{currently not used.
    }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  For each descriptor, a cumulative probit model is fitted and the
  relevant information extracted and from the fits.

  The function has \code{print} and \code{plot} methods.
}
\value{
  An object of class \code{opair} with components:
  \item{coefficients}{a table (\code{data.frame}) of d-prime estimates,
    standard errors, confidence limits and p-values for difference and
    equivalence tests.}
  \item{convergence}{a vector indicating if the by-descriptor cumulative
    link models converged: zero indicates convergence.}
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
  %%  ~~who you are~~
  Rune Haubo B Christensen
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{plot.opair}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Use the NV data set:
names(NV)
## get vector of descriptor/attribute names:
(attr.names <- names(NV)[3:ncol(NV)])
## Fit Thurstonian ordinal paired comparison model:
opair(descriptors=NV[attr.names],
      products=NV$Samples, d.equiv=1)

## Use the BG data set:
names(BG)
## get vector of descriptor/attribute names:
(attr.names <- names(BG)[3:ncol(BG)])
## Fit Thurstonian ordinal paired comparison model:
opair(descriptors=BG[attr.names],
      products=BG$Samples, d.equiv=1, conf.level=.90)

## Abbreviate names:
opair(BG[attr.names], BG$Samples, d.equiv=1,
conf.level=.90, abbreviate=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
