normalPwr <-
  function(pdA, pd0 = 0, sample.size, alpha = 0.05, pGuess = 1/2,
           test = c("difference", "similarity"))
{
  test <- match.arg(test)
  stopifnot(is.numeric(pdA) && length(pdA) == 1 &&
            pdA >= 0 && pdA <= 1)
  stopifnot(is.numeric(pd0) && length(pd0) == 1 &&
            pd0 >= 0 && pd0 <= 1)
  stopifnot(is.numeric(sample.size) && length(sample.size) == 1 &&
            isTRUE(all.equal(round(sample.size), sample.size)) &&
            sample.size > 0)
  sample.size <- as.integer(round(sample.size))
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(is.numeric(pGuess) && length(pGuess) == 1 &&
            pGuess >= 0 && pGuess < 1)
  if(test == "difference" && pdA <= pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA >= pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  pc0 <- pd2pc(pd=pd0, Pguess=pGuess)
  pcA <- pd2pc(pd=pdA, Pguess=pGuess)
  sigma0 <- sqrt(pc0*(1 - pc0)/sample.size)
  sigmaA <- sqrt(pcA*(1 - pcA)/sample.size)
  if(test == "difference") {
    lambda <- (qnorm(1 - alpha) * sigma0 + pc0 - pcA) / sigmaA
    pwr <- pnorm(lambda, lower.tail = FALSE)
  }
  else if(test == "similarity") {
    lambda <- (qnorm(alpha) * sigma0 + pc0 - pcA) / sigmaA
    pwr <- pnorm(lambda, lower.tail = TRUE)
  }
  else
    stop("'test' not recognized")
  return(pwr)
}

discrimPwr <-
  function(pdA, pd0 = 0, sample.size, alpha = 0.05, pGuess = 1/2,
           test = c("difference", "similarity"),
           statistic = c("exact", "normal"))
{
  ## match and test arguments:
  test <- match.arg(test)
  stat <- match.arg(statistic)
  ss <- sample.size
  stopifnot(is.numeric(pdA) && length(pdA) == 1 &&
            pdA >= 0 && pdA <= 1)
  stopifnot(is.numeric(pd0) && length(pd0) == 1 &&
            pd0 >= 0 && pd0 <= 1)
  stopifnot(is.numeric(sample.size) && length(sample.size) == 1 &&
            isTRUE(all.equal(round(sample.size), sample.size)) &&
            sample.size > 0)
  sample.size <- as.integer(round(sample.size))
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(is.numeric(pGuess) && length(pGuess) == 1 &&
            pGuess >= 0 && pGuess < 1)

  ## Get pc from pdA and pGuess:
  pc <- pd2pc(pdA, pGuess)
  if(stat == "normal") {
    pwr <- normalPwr(pdA = pdA, pd0 = pd0, sample.size = ss,
                     alpha = alpha, pGuess = pGuess, test = test)
    return(pwr)
  } ## stat == "exact":
  ## critical value in one-tailed binomial test:
  xcr <- findcr(sample.size=ss, alpha=alpha, p0=pGuess, pd0=pd0,
                test=test) 
  ## compute power of the test from critical value:
  if(test == "difference") {
    xcr <- delimit(xcr, lower = 1, upper = ss + 1)
    power <- 1 - pbinom(q = xcr - 1, size = ss, prob = pc)
  }
  else if(test == "similarity") {
    xcr <- delimit(xcr, lower = 0, upper = ss)
    power <- pbinom(q = xcr, size = ss, prob = pc)
  }
  else ## should never happen
    stop("'test' not recognized")
  return(power)
}

d.primePwr <-
  function(d.primeA, d.prime0 = 0, sample.size, alpha = 0.05,
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle"), 
           test = c("difference", "similarity"),
           statistic = c("exact", "normal"))
{
  ## Convenience function that simply modifies some arguments and
  ## calls discrimPwr
  newCall <- call <- match.call()
  method <- match.arg(method)
  stopifnot(length(d.primeA) == 1 && is.numeric(d.primeA) &&
            d.primeA >= 0)
  stopifnot(length(d.prime0) == 1 && is.numeric(d.prime0) &&
            d.prime0 >= 0)
  pdA <- coef(rescale(d.prime = d.primeA, method = method))$pd
  pd0 <- coef(rescale(d.prime = d.prime0, method = method))$pd
  newCall$method <- newCall$d.primeA <- newCall$d.prime0 <- NULL
  newCall$pGuess <-
    ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  newCall$pdA <- pdA
  newCall$pd0 <- pd0
  newCall[[1]] <- as.name("discrimPwr")
  return(eval.parent(newCall))
}

normalSS <-
  function(pdA, pd0 = 0, target.power = 0.9, alpha = 0.05,
           pGuess = 1/2, test = c("difference", "similarity"))
{
  test <- match.arg(test)
  stopifnot(pdA >= 0 && pdA <= 1)
  stopifnot(pd0 >= 0 && pdA <= 1)
  stopifnot(target.power > 0 && target.power < 1)
  stopifnot(alpha > 0 && alpha < 1)
  stopifnot(pGuess >= 0 && pGuess < 1)
  if(test == "difference" && pdA <= pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA >= pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  pc0 <- pd2pc(pd0, pGuess)
  pcA <- pd2pc(pdA, pGuess)
  s0 <- sqrt(pc0*(1 - pc0))
  sA <- sqrt(pcA*(1 - pcA))
  ## approximate sample size:
  if(test == "difference") {
    stopifnot(pdA > pd0)
    n <- ((sA * qnorm(1 - target.power) - s0 * qnorm(1 - alpha)) /
          (pc0 - pcA))^2
  }
  else if(test == "similarity") {
    stopifnot(pdA < pd0)
    n <- ((sA * qnorm(target.power) - s0 * qnorm(alpha)) /
          (pc0 - pcA))^2
  }
  return(ceiling(n))
}

discrimSS <-
  function(pdA, pd0 = 0, target.power = 0.90, alpha = 0.05,
           pGuess = 1/2, test = c("difference", "similarity"),
           statistic = c("exact", "normal")) 
{
  test <- match.arg(test)
  call <- match.call()
  stat <- match.arg(statistic)
  stopifnot(length(pdA) == 1 && length(pd0) == 1 &&
            length(target.power) == 1 && length(alpha) == 1 && 
            length(pGuess) == 1)
  if(pdA < 0 | pdA > 1)
    stop("'pdA' has to be between zero and one")
  if(alpha <= 0 | alpha >= 1)
    stop("'alpha' has to be between zero and one")
  if(target.power <= 0 | target.power >= 1)
    stop("'target.power' has to be between zero and one")
  if(pd0 < 0 | pd0 > 1)
    stop("'pd0' has to be between zero and one")
  if(test == "difference" && pdA <= pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA >= pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  ssN <- normalSS(pdA = pdA, pd0 = pd0, target.power = target.power,
                  alpha = alpha, pGuess = pGuess, test = test)
  ssN <- ssN - 1
  if(stat == "normal")
    return(ssN + 1)
  ## if(stat == "exact"):
  if(ssN <= 50)
    ssTry <- 1
  if(ssN > 1e4) {
    warning(paste("sample size probably > 1e4 and 'exact' option is",
                  "not available\n",
                  "using normal approximation instead"))
    return(ssN)
  }
  if(ssN > 50 && ssN <= 1e4) {
    for(ssTry in c(floor(ssN * (9:1 / 10)), 1) ) {
      tryPwr <-
        discrimPwr(pdA = pdA, pd0 = pd0, sample.size = ssTry,
                   alpha = alpha, pGuess = pGuess, test = test)
      if(tryPwr < target.power)
        break
    }
    if(tryPwr > target) stop("Failed to find admissible sample size")
  }
  while (discrimPwr(pdA = pdA, pd0 = pd0, sample.size = ssTry, 
                    alpha = alpha, pGuess = pGuess, 
                    test = test) < target.power)
    ssTry <- ssTry + 1
  return(ssTry)
}

d.primeSS <- 
  function(d.primeA, d.prime0 = 0, target.power = 0.90, alpha = 0.05,
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle"), 
           test = c("difference", "similarity"),
           statistic = c("exact", "normal")) 
{
  ## Convenience function that simply modifies some arguments and
  ## calls discrimSS
  newCall <- call <- match.call()
  method <- match.arg(method)
  stopifnot(length(d.primeA) == 1 && is.numeric(d.primeA) &&
            d.primeA >= 0)
  stopifnot(length(d.prime0) == 1 && is.numeric(d.prime0) &&
            d.prime0 >= 0)
  pdA <- coef(rescale(d.prime = d.primeA, method = method))$pd
  pd0 <- coef(rescale(d.prime = d.prime0, method = method))$pd
  newCall$method <- newCall$d.primeA <- newCall$d.prime0 <- NULL
  newCall$pGuess <-
    ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  newCall$pdA <- pdA
  newCall$pd0 <- pd0
  newCall[[1]] <- as.name("discrimSS")
  return(eval.parent(newCall))
}

