## Various functions used to calculate bias of d-prime in sensory
## discrimination protocols.

pbelow <-
    function(d.prime, size,
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"))
### Probability that an observation will fall below pG:
### Prob(X <= floor(n*pG))
{
    method <- match.arg(method)
    pA <- psyfun(d.prime, method)
    pG <- if(method %in% c("duotrio", "twoAFC")) 1/2 else 1/3
    xsmall <- floor(size * pG)
    pbinom(xsmall, size, pA)
}

plast <-
    function(d.prime, size,
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"))
### Probability that an observation will fall in x=size
{
    pA <- psyfun(d.prime, method)
    dbinom(size, size, prob=pA)
}

dbias <-
    function(d.prime, size,
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
             last=c("add.last", "ignore"))
### Compute bias in d.prime.
### Note that X=size => d.prime=Inf, which is handled in two ways:
### 1) Adding Prob(X = size) to Prob(X = size - 1), or
### 2) Ignoring Prob(X = size)
{
    method <- match.arg(method)
    last <- match.arg(last)
    pA <- psyfun(d.prime, method)
    pG <- if(method %in% c("duotrio", "twoAFC")) 1/2 else 1/3
    X <- pdens(pA=pA, n=size, pg=pG)
    m <- nrow(X)
    dp <- psyinv(X[-m, "phat"], method)
    if(last == "add.last") {
        dens <- X[, "dens"]
        dens <- c(dens[seq_len(m-2)], sum(rev(dens)[1:2]))
    } else if(last == "ignore") {
        dens <- X[-m, "dens"]
    } else stop("'last' not recognized")
    bias <- sum(dp * dens) - d.prime
    bias
}

## Could we also evaluate variance and MSE of the estimator?

pbias <- function(pA, n, pg=1/2)
### Compute bias in binomial p.
{
    x <- 0:n
    phat <- delimit(x=x/n, lower=pg, upper=1)
    d <- dbinom(x=x, size=n, prob=pA)
    bias <- sum(phat * d) - pA
    bias
}

pdens <- function(pA, n, pg=1/2)
### Compute density of outcomes.
{
    xsmall <- floor(n*pg)
    dens <- c(pbinom(xsmall, n, pA),
              dbinom(x=(xsmall+1):n, size=n, prob=pA))
    stopifnot(abs(sum(dens) - 1) < sqrt(.Machine$double.eps))
    cbind(x=xsmall:n, phat=(xsmall:n)/n, dens=dens)
}

