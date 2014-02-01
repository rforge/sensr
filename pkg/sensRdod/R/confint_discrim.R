expCIwidthd <-
    function(d.prime, n, level=0.95,
             method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
             statistic=c("likelihood", "exact", "score", "Waldp", "Waldd"), ...)
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
                   statistic=c("likelihood", "exact", "score", "Waldp"))
### Compute the CI for the binomial p.
{
    stat <- match.arg(statistic)
    if(stat == "likelihood") {
        ci <- lrci(x=x, n=n, level=level)
    } else if(stat == "exact") {
        ci <- binom.test(x=x, n=n, conf.level=level)$conf.int
    } else if(stat == "Waldp") {
        ci <- waldcip(x=x, n=n, level=level)
    } else if(stat == "score") {
        ci <- prop.test(x=x, n=n, conf.level=level)$conf.int
    }
    ci <- setNames(c(ci), c("lower", "upper"))
    attr(ci, "level") <- level
    ci
}

getCId <- function(x, n, level=.95, plow=c("pGuess", "NA"),
                   method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
                   statistic=c("likelihood", "exact", "score", "Waldp", "Waldd"))
### Compute the CI for d-prime.
{
    method <- match.arg(method)
    stat <- match.arg(statistic)
    if(stat == "Waldd") {
        ci <- walddci(x=x, n=n, level=level, plow=plow, method=method)
        ci <- delimit(ci, lower=0)
    } else {
        ci <- getCIp(x=x, n=n, level=level, statistic=stat)
        ci <- psyinv(delimit(ci, lower=0, upper=1), method=method)
    }
    ci <- setNames(c(ci), c("lower", "upper"))
    attr(ci, "level") <- level
    ci
}

CIdistp <-
    function(pA, n, pG=1/2, level=0.95,
             statistic=c("likelihood", "exact", "score", "Waldp"), ...)
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
             statistic=c("likelihood", "exact", "score", "Waldp", "Waldd"), ...)
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
    attr(ci, "level") <- level
    ci
}

waldcip <- function(x, n, level=0.95) {
    p <- x/n
    se <- sqrt(p * (1 - p) / n)
    waldci(par=p, se=se, level=level)
}

lrci <-
    function(x, n, neval=100, level=.95)
{
    drop(confint(profBinom(x, n, nProf=neval), level=level))
}

lrci_d <-
    function(x, n, method=c("duotrio", "tetrad", "threeAFC", "twoAFC", "triangle"),
             neval=100, level=.95)
{
    method <- match.arg(method)
    ci <- lrci(x=x, n=n, neval=neval, level=level)
    ci[] <- psyinv(ci, method=method)
    ci
}
