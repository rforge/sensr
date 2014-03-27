expCIwidthd <-
    function(d.prime, n, level=0.95,
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
             statistic=c("likelihood", "exact", "score", "waldp", "waldd"), ...)
### Compute the expected/average and most probable CI width for d-prime.
{
    call <- match.call()
    call[[1]] <- as.name("CIdistd")
    XX <- eval.parent(call)
    XX <- cbind(XX, width=apply(XX[, c("lower", "upper")], 1, diff))
    XX <- XX[is.finite(XX[, "width"]), ]
    ## Expected and mostProb. among the finite length intervals:
    expCIw <- sum(apply(XX[, c("dens", "width")], 1, prod))
    mpCIw <- XX[which.max(XX[, "dens"]), "width"]
    c("expected"=as.vector(expCIw), "most.prob"=as.vector(mpCIw))
}

getCIp <- function(x, n, level=.95,
                   statistic=c("likelihood", "exact", "score",
                   "waldp", "blaker", "AC"))
### Compute the CI for the binomial p.
{
    stat <- match.arg(statistic)
    if(stat == "likelihood") {
        ci <- lici(x=x, n=n, level=level)
    } else if(stat == "exact") {
        ci <- binom.test(x=x, n=n, conf.level=level)$conf.int
    } else if(stat == "waldp") {
        ci <- waldcip(x=x, n=n, level=level)
    } else if(stat == "score") {
        ci <- prop.test(x=x, n=n, conf.level=level)$conf.int
    } else if(stat == "AC") {
        ci <- ACci(x=x, n=n, level=level)
    } else if(stat == "blaker") {
        ci <- blakerci(x=x, n=n, level=level)
    }
    ci <- setNames(c(ci), c("lower", "upper"))
    ## attr(ci, "level") <- level
    ci
}

getCId <- function(x, n, level=.95, plow=c("pGuess", "NA"),
                   method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
                   statistic=c("likelihood", "exact", "score",
                   "blaker", "waldp", "waldd"))
### Compute the CI for d-prime.
{
    method <- match.arg(method)
    stat <- match.arg(statistic)
    if(stat == "waldd") {
        ci <- walddci(x=x, n=n, level=level, plow=plow, method=method)
        ci <- delimit(ci, lower=0)
    } else {
        ci <- getCIp(x=x, n=n, level=level, statistic=stat)
        ci <- psyinv(delimit(ci, lower=0, upper=1), method=method)
    }
    ci <- setNames(c(ci), c("lower", "upper"))
    ## attr(ci, "level") <- level
    ci
}

CIdistp <-
    function(pA, n, pG=1/2, level=0.95,
             statistic=c("likelihood", "exact", "score", "waldp"), ...)
### Computes the distribution of CIs for a binomial p.
{
    stat <- match.arg(statistic)
    XX <- pdens(pA, n=n, pg=pG)
    CI <- do.call(rbind, lapply(XX[, "x"], getCIp, n=n, statistic=stat))
    colnames(CI) <- c("lower", "upper")
    cbind(XX, CI)
}

CIdistd <-
    function(d.prime, n, level=0.95,
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
             statistic=c("likelihood", "exact", "score", "waldp", "waldd"), ...)
### Computes the distribution of CIs for d-prime.
{
    method <- match.arg(method)
    stat <- match.arg(statistic)
    pA <- psyfun(d.prime, method)
    pG <- if(method %in% c("duotrio", "twoAFC")) 1/2 else 1/3
    XX <- pdens(pA, n=n, pg=pG)
    CI <- lapply(XX[, "x"], function(x) {
        getCId(x=x, n=n, level=level, plow="pGuess",
               method=method, statistic=statistic)
    })
    CI <- do.call(rbind, CI)
    colnames(CI) <- c("lower", "upper")
    cbind(XX, CI)
}

walddci <-
    function(x, n, level=0.95, plow=c("pGuess", "NA"),
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"))
{
    method <- match.arg(method)
    plow <- match.arg(plow)
    a <- 1 - level
    pG <- if(method %in% c("duotrio", "twoAFC")) 1/2 else 1/3
    phat <- if(plow == "pGuess") delimit(x/n, lower=pG) else x/n
    se.p <- sqrt(phat * (1 - phat) / n)
    dp <- psyinv(phat, method=method)
    se.dp <- se.p / psyderiv(dp, method = method)
    ci <- waldci(dp, se.dp, level=level)
    ci[] <- delimit(ci, lower=0)
    ci
}

waldci <- function(par, se, level=0.95) {
    a <- 1-level
    ci <- par + c(-1, 1) * qnorm(1 - a/2) * se
    ci <- setNames(ci, c("lower", "upper"))
    ## attr(ci, "level") <- level
    ci
}

waldcip <- function(x, n, level=0.95) {
    p <- x/n
    se <- sqrt(p * (1 - p) / n)
    waldci(par=p, se=se, level=level)
}

## lrci <-
##     function(x, n, neval=100, level=.95)
## {
##     drop(confint(profBinom(x, n, nProf=neval), level=level))
## }

lici_d <-
    function(x, n, method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
             level=.95)
{
    method <- match.arg(method)
    ## ci <- lrci(x=x, n=n, neval=neval, level=level)
    ci <- lici(x=x, n=n, level=level)
    ci[] <- psyinv(ci, method=method)
    ci
}

lici <- function(x, n, level=0.95) {
    rellike <- function(theta, cutoff=0)
        dbinom(x, n, theta) / ML - cutoff
    co <- exp(-qchisq(level, df=1)/2) ## cutoff
    phat <- x/n
    ML <- dbinom(x, n, phat)
    lwr <- 0
    upr <- 1
    if(phat > 0)
        lwr <- uniroot(rellike, c(lwr, phat), cutoff=co)$root
    if(phat < 1)
        upr <- uniroot(rellike, c(phat, upr), cutoff=co)$root
    c("lwr"=lwr, "upr"=upr)
}

## waldci <- function(x, n, level=0.95) {
##     phat <- x/n
##     se.phat <- sqrt(phat * (1 - phat) / n)
##     ci <- phat + c(-1, 1) * qnorm(1-alpha/2) * se.phat
##     ci[ci < 0] <- 0
##     ci[ci > 1] <- 1
##     ci
## }

ACci <- function(x, n, level=0.95) {
    waldcip(x+2, n+4, level=level)
}

scoreci <- function(x, n, level=0.95) {
    ci <- suppressWarnings(
        c(prop.test(x=x, n=n, alternative="two.sided",
                    conf.level=level, correct=FALSE)$conf.int)
        )
    ci
}

exactci <- function(x, n, level=0.95) {
    ci <- c(binom.test(x=x, n=n, alternative="two.sided",
                 conf.level=level)$conf.int)
}

blakerci <- function(x, n, level=0.95) {
    acceptbin <- function(x, n, p) {
        ## Computes the acceptability of p when x is observed and X is
        ## Bin(n,p)
        p1 <- 1 - pbinom(x - 1, n, p)
        p2 <- pbinom(x, n, p)
        a1 <- p1 + pbinom(qbinom(p1, n, p) - 1, n, p)
        a2 <- p2 + 1 - pbinom(qbinom(1 - p2, n, p), n, p)
        min(a1,a2)
    }
    acc <- function(theta, alpha)
        acceptbin(x, n, p=theta) - alpha
    alpha <- 1 - level
    phat <- x/n
    lwr <- 0
    upr <- 1
    if(phat > 0)
        lwr <- uniroot(acc, c(lwr, phat), alpha=alpha)$root
    if(phat < 1)
        upr <- uniroot(acc, c(phat, upr), alpha=alpha)$root
    c("lwr"=lwr, "upr"=upr)
}

binomci <-
    function(x, n, level=0.95,
             statistic=c("likelihood", "waldp", "AC", "score", "exact", "blaker"))
{
    stopifnot(length(n) == 1L)
    stopifnot(length(x) >= 1L)
    stopifnot(all(x >= 0) && all(x <= n) &&
              all(abs(round(x) - x) < 1e-4))
    x <- as.integer(round(x))
    statistic <- match.arg(statistic)
    cifun <- switch(statistic,
                    "likelihood"=lici,
                    "waldp"=waldcip,
                    "AC"=ACci,
                    "score"=scoreci,
                    "exact"=exactci,
                    "blaker"=blakerci)
    ci.list <- lapply(x, cifun, n=n, level=level)
    ci.mat <- do.call(rbind, ci.list)
    dimnames(ci.mat) <- list(as.character(x), c("lwr", "upr"))
    ci.mat
}

binomci.d <-
    function(x, n, level=0.95, plow=c("pGuess", "NA"),
             statistic=c("likelihood", "waldp", "AC", "score", "exact",
             "blaker", "waldd"),
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle"))
{
    stopifnot(length(n) == 1L)
    stopifnot(length(x) >= 1L)
    stopifnot(all(x >= 0) && all(x <= n) &&
              all(abs(round(x) - x) < 1e-4))
    x <- as.integer(round(x))
    statistic <- match.arg(statistic)
    method <- match.arg(method)
    plow <- match.arg(plow)
    if(statistic == "waldd") {
        ci.list <- lapply(x, walddci, n=n, level=level, plow=plow, method=method)
        ci.mat <- do.call(rbind, ci.list)
        ci.mat[] <- delimit(c(ci.mat), lower=0)
    } else {
        ## cifun <- switch(statistic,
        ##                 "like"=lici,
        ##                 "wald"=waldcip,
        ##                 "AC"=ACci,
        ##                 "score"=scoreci,
        ##                 "exact"=exactci,
        ##                 "blaker"=blakerci)
        ## ci.list <- lapply(x, cifun, n=n, level=level)
        ## ci.mat <- do.call(rbind, ci.list)
        ci.mat <- binomci(x=x, n=n, level=level, statistic=statistic)
        ci.mat[] <- psyinv(delimit(c(ci.mat), lower=0, upper=1),
                           method=method)
    }
    dimnames(ci.mat) <- list(as.character(x), c("lwr", "upr"))
    ci.mat
}


covProb <-
    function(theta, n, level=0.95,
             method=c("likelihood", "waldp", "AC", "score", "exact", "blaker"))
{
    ## Vectorized over theta.
    stopifnot(length(n) == 1L)
    stopifnot(abs(round(n) - n) < 1e-4)
    n <- as.integer(round(n))
    stopifnot(length(theta) >= 1L)
    stopifnot(all(theta <=1 & theta >= 0))
    method <- match.arg(method)
    xvec <- 0:n
    ## dens <- dbinom(xvec, n, theta)
    cifun <- switch(method,
                    "likelihood"=lici,
                    "waldp"=waldcip,
                    "AC"=ACci,
                    "score"=scoreci,
                    "exact"=exactci,
                    "blaker"=blakerci)
    CI <- t(sapply(xvec, cifun, n=n, level=level))
    cp <- sapply(theta, function(th) {
        ok <- CI[, 1] <= th & th <= CI[, 2]
        sum(dbinom(xvec, n, th)[ok]) ## Coverage probability
    })
    ## ok <- CI[, 1] <= theta & theta <= CI[, 2]
    ## cp <- sum(dens[ok]) ## Coverage probability
    cp
}

covProb(.5, 10, method="like")
covProb(.6, 10, method="like")

covProb(c(.5, .6), 10, method="like")

ciWidth <-
    function(theta, n, level=0.95, type=c("expected", "median"),
             method=c("likelihood", "waldp", "AC", "score", "exact", "blaker"))
{
    ## Vectorized over theta.
    stopifnot(length(n) == 1L)
    stopifnot(abs(round(n) - n) < 1e-4)
    stopifnot(length(theta) >= 1L)
    stopifnot(all(theta <=1 & theta >= 0))
    method <- match.arg(method)
    type <- match.arg(type)
    ## if(type == "median") stop("'type = median' is not yet implemented")
    n <- as.integer(round(n))
    xvec <- 0:n
    ## dens <- dbinom(xvec, n, theta)
    cifun <- switch(method,
                    "likelihood"=lici,
                    "waldp"=waldcip,
                    "AC"=ACci,
                    "score"=scoreci,
                    "exact"=exactci,
                    "blaker"=blakerci)
    CI <- t(sapply(xvec, cifun, n=n, level=level))
    CI.width <- abs(apply(CI, 1, diff))
    if(type == "expected") {
        ll <- vapply(theta, function(th) {
            weighted.mean(CI.width, dbinom(xvec, n, prob=th))
        }, numeric(1))
    } else { # type == "median"
        ll <- vapply(theta, function(th) {
            CI.width[which.max(dbinom(xvec, n, prob=th))[1L]]
        }, numeric(1))
    }
    ll
}

ciWidth.d <-
    function(d.prime, n, level=0.95, type=c("expected", "median"),
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
             statistic=c("likelihood", "waldp", "score",
             "exact", "blaker", "waldd"))
{
    ## Vectorized over d.prime.
    stopifnot(length(n) == 1L)
    stopifnot(abs(round(n) - n) < 1e-4)
    stopifnot(length(d.prime) >= 1L)
    stopifnot(all(d.prime >= 0))
    method <- match.arg(method)
    type <- match.arg(type)
    stat <- match.arg(statistic)
    n <- as.integer(round(n))
    xvec <- 0:n
    theta <- psyfun(d.prime, method=method)
    ## Compute all CI and their widths:
    CI <- t(sapply(xvec, getCId, n=n, level=level, method=method,
                   statistic=stat))
    keep <- apply(CI, 1, function(x) all(is.finite(x)))
    xvec <- xvec[keep]
    CI.width <- abs(apply(CI[keep, ], 1, diff))
    ## Average CI-width
    if(type == "expected") {
        ll <- vapply(theta, function(th) {
            weighted.mean(CI.width, dbinom(xvec, n, prob=th))
        }, numeric(1))
        ## Most probable CI-width:
    } else { # type == "median"
        ll <- vapply(theta, function(th) {
            CI.width[which.max(dbinom(xvec, n, prob=th))[1L]]
        }, numeric(1))
    }
    ll
}

covProb.waldd <-
    function(d.prime, n, level=0.95, plow=c("pGuess", "NA"),
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"))
{
    ## Vectorized over theta.
    stopifnot(length(n) == 1L)
    stopifnot(abs(round(n) - n) < 1e-4)
    n <- as.integer(round(n))
    stopifnot(length(d.prime) >= 1L)
    stopifnot(all(d.prime >= 0))
    method <- match.arg(method)
    theta <- psyfun(d.prime, method=method)
    xvec <- 0:n
    CI <- do.call(rbind, lapply(xvec, walddci, n=n, level=level,
                                plow=plow, method=method))
    keep <- complete.cases(CI)
    xvec.keep <- xvec[keep]
    CI.keep <- CI[keep, , drop=FALSE]
    norm.const <- sum(dens.keep)
    cp <- sapply(seq_along(d.prime), function(i) {
        ok <- CI.keep[, 1] <= d.prime[i] & d.prime[i] <= CI.keep[, 2]
        dens <- dbinom(xvec.keep, n, theta[i])
        sum(dens[ok]) ## Coverage probability
    })
    cp
}

